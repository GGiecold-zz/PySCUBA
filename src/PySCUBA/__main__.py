#!/usr/bin/env python


# PySCUBA/src/PySCUBA/__main__.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from os import getcwd, path, remove
import Queue
import sys

try:
    import igraph
except ImportError, e:
    pass

import numpy as np
from PIL import Image, ImageQt
from PyQt4 import QtCore, QtGui
from sklearn.preprocessing import StandardScaler
import wand.image

from .Gap_stats import gap_stats
from .Preprocessing import cytometry_preprocess, PCR_preprocess, RNASeq_preprocess
from . import PySCUBA_design, SCUBA_core


def plot_tree(cluster_indices, parent_clusters, output_directory = None):
    """Display a bifurcation tree.
    """

    if igraph not in sys.modules:
        return

    if output_directory is None:
        output_directory = getcwd()
        
    vertex_sizes = np.bincount(cluster_indices)
    N_vertices = vertex_sizes.size
    vertex_sizes = np.divide(vertex_sizes, float(np.sum(vertex_sizes)))
    vertex_sizes *= 100 * N_vertices
    vertex_sizes += 40 + (N_vertices / 3)

    tree = igraph.Graph()
    tree.add_vertices(N_vertices)
    
    cluster_tally = 0
    for k, v in parent_clusters.items():
        if k > 0:
            tree.add_edges(zip(v, xrange(cluster_tally, cluster_tally + len(v))))
        cluster_tally += len(v)
        
    tree.vs['label'] = xrange(N_vertices)
    layout = tree.layout('fr')
    name = path.join(output_directory, 'SCUBA_tree.pdf')
    igraph.plot(tree, name, bbox = (200 * N_vertices, 200 * N_vertices), margin = 250,
                layout = layout, edge_width = [7] * (N_vertices - 1), 
                vertex_label_dist = 0, vertex_label_size = 30, 
                vertex_size = vertex_sizes.tolist())
    

def one_to_max(array_in):
    """Alter a vector of cluster labels to a dense mapping. 
        Given that this function is herein always called after passing 
        a vector to the function checkcl, one_to_max relies on the assumption 
        that cluster_run does not contain any NaN entries.

    Parameters
    ----------
    array_in : a list or one-dimensional array
        The list of cluster IDs to be processed.
    
    Returns
    -------
    result : one-dimensional array
        A massaged version of the input vector of cluster identities.
    """
    
    x = np.asanyarray(array_in)
    N_in = x.size
    array_in = x.reshape(N_in)    

    sorted_array = np.sort(array_in)
    sorting_indices = np.argsort(array_in)

    last = np.nan
    current_index = -1
    for i in xrange(N_in):
        if last != sorted_array[i] or np.isnan(last):
            last = sorted_array[i]
            current_index += 1

        sorted_array[i] = current_index

    result = np.empty(N_in, dtype = int)
    result[sorting_indices] = sorted_array

    return result
       

class WorkerThread(QtCore.QThread):

    def __init__(self, result_queue, data_type, data_path, cluster_mode, log_mode,
                 pseudotime_mode, pcv_method, anchor_gene,
                 exclude_marker_names):
        super(WorkerThread, self).__init__()
        
        self.result_queue = result_queue
        self.data_type = str(data_type)
        self.data_path = data_path
        
        cluster_mode = str(cluster_mode).lower()
        self.cluster_mode = None if (cluster_mode == 'none') else cluster_mode
            
        self.log_mode = log_mode
        self.pseudotime_mode = pseudotime_mode
        self.pcv_method = pcv_method
        self.anchor_gene = anchor_gene
        self.exclude_marker_names = exclude_marker_names
        
    def __del__(self):
        self.wait()
        
    def run(self):
        preprocessing_fcts = [cytometry_preprocess, PCR_preprocess,
                              RNASeq_preprocess]
        data_type_dict = {'cytometry': 0, 'PCR': 1, 'RNASeq': 2}
        cell_IDs, data, markers, cell_stages, data_tag, \
        output_directory = preprocessing_fcts[data_type_dict[self.data_type]](
            self.data_path, self.log_mode, self.pseudotime_mode,
            self.pcv_method, self.anchor_gene, self.exclude_marker_names)
    
        cell_stages = 1 + one_to_max(cell_stages)
        data = StandardScaler(with_std = False).fit_transform(data)
    
        if self.cluster_mode in {'pca', 'pca2'}:    
            PCA_components, data = SCUBA_core.PCA_analysis(data, self.cluster_mode,
                cell_stages if (self.cluster_mode == 'pca2') else None)
    
        centroid_coords, cluster_indices, \
        parent_clusters = SCUBA_core.initialize_tree(data, cell_stages)
        centroid_coords, cluster_indices, \
        parent_clusters = SCUBA_core.refine_tree(data, centroid_coords,
            cluster_indices, parent_clusters, cell_stages, output_directory)
    
        plot_tree(cluster_indices, parent_clusters, output_directory)
        
        if self.cluster_mode in {'pca', 'pca2'}:
            weights = PCA_components
        else:
            weights = None
    
        bifurcation_info, bifurcation_axes, \
        bifurcation_projections = SCUBA_core.bifurcation_direction(data, cell_IDs,
            markers, parent_clusters, centroid_coords, output_directory,
            weights)
                
        if bifurcation_info:
            data_per_split, parameters_per_split = SCUBA_core.bifurcation_analysis(
                cluster_indices, bifurcation_info, bifurcation_axes,
                bifurcation_projections, output_directory)
                
        self.result_queue.put(output_directory)
        self.result_queue.task_done()
        
        return


class LoadImageThread(QtCore.QThread):

    def __init__(self, source_file):
        super(LoadImageThread, self).__init__()
        self.source_file = source_file
        
    def __del__(self):
        self.wait()
        
    def run(self):
        self.emit(QtCore.SIGNAL("showImage(QString)"), self.source_file)


class PySCUBApp(QtGui.QMainWindow, PySCUBA_design.Ui_MainWindow):

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)
        
        self.cancelButton.setEnabled(False)
        self.okButton.setEnabled(False)

        self.data_path = './'
        self.selectDatasetButton.clicked.connect(self.selectDataset)
        
        self.log_mode = True
        self.logCheckBox.stateChanged.connect(self.logStateChanged)
        
        self.pseudotime_mode = True
        self.pseudotimeCheckBox.stateChanged.connect(
            self.pseudotimeStateChanged)
        
        self.pcv_method = 'Rprincurve' 
        self.anchor_gene = None 
        self.exclude_marker_names = None
        
        self.result_queue = Queue.Queue()
        self.okButton.clicked.connect(self.buttonClicked)
        self.okButton.clicked.connect(self.OK)
        self.cancelButton.clicked.connect(self.buttonClicked)
        
        self.zoom = 0
        self.pixMap = QtGui.QPixmap()
        self.displayFileButton.setEnabled(False)
        self.displayFileButton.clicked.connect(self.selectDisplay)

    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(self, 'Message', 
            "Are you sure to quit?", QtGui.QMessageBox.Yes 
            | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.close()

    def selectDataset(self):
        dataFileDialog = QtGui.QFileDialog(self)
        self.data_path = str(dataFileDialog.getOpenFileName())
        self.statusbar.showMessage("{0} ready to be "
            "analyzed".format(path.basename(self.data_path)))
        self.cancelButton.setEnabled(True)
        self.okButton.setEnabled(True)
    
    def logStateChanged(self, int):
        if self.logCheckBox.isChecked():
            self.log_mode = True
        else:
            self.log_mode = False
    
    def pseudotimeStateChanged(self, int):
        if self.pseudotimeCheckBox.isChecked():
            self.pseudotime_mode = True
        else:
            self.pseudotime_mode = False
    
    def buttonClicked(self):
        sender = self.sender()
        self.statusbar.showMessage(sender.text() + " was pressed.")
        self.button_clicked = sender.text()
        
    def OK(self):
        self.statusbar.showMessage('Work in progress...')
        
        self.worker_thread = WorkerThread(self.dataTypeComboBox.currentText(),
            self.data_path, self.clusterModeComboBox.currentText(),
            self.log_mode, self.pseudotime_mode, self.pcv_method,
            self.anchor_gene, self.exclude_marker_names)
        self.connect(self.worker_thread, QtCore.SIGNAL("update(QString)"),
            self.worker_thread.run)
        self.connect(self.worker_thread, QtCore.SIGNAL("finished()"), self.doneRunning)
        self.worker_thread.start()
        self.cancelButton.setEnabled(True)
        self.okButton.setEnabled(False)
    
    def doneRunning(self):
        if self.button_clicked == 'Cancel':
            self.cancelRunning()
        else:
            self.cancelButton.setEnabled(False)
            self.okButton.setEnabled(False)
            self.displayFileButton.setEnabled(True)
            self.directory = self.result_queue.get()
            self.statusbar.showMessage("PySCUBA has completed the "
                "analysis of your data.")
            QtGui.QMessageBox.information(self, "Status Message", 
                "Mission accomplished!")
            
    def cancelRunning(self):
        self.cancelButton.setEnabled(False)
        self.okButton.setEnabled(False)
        self.worker_thread.terminate()
        self.statusbar.showMessage("PySCUBA was interrupted!")
        QtGui.QMessageBox.information(self, "Status Message", 
            "PySCUBA was interrupted!")
    
    def selectDisplay(self):
        filters = 'Images (*.jpg *.pdf *.png)'
        select_filters = 'Images (*.jpg *.pdf *.png)'
        source_file = QtGui.QFileDialog.getOpenFileName(self, 
            'Select file to display', self.directory, filters, select_filters)
        
        self.load_image_thread = LoadImageThread(source_file)
        self.connect(self.load_image_thread, QtCore.SIGNAL("showImage(QString)"),
            self.showImage)
        self.load_image_thread.start()
        
    def zoomFactor(self):
        return self.zoom
        
    def wheelEvent(self, event):
        if not self.pixMap.isNull():
            if event.delta() < 0:
                factor = 0.8
                self.zoom -= 1
            else:
                factor = 1.25
                self.zoom += 1
            if self.zoom < 0:
                self.zoom = 0
            elif self.zoom == 0:
                self.fitInView()
            else:
                self.graphicsView.scale(factor, factor)
        else:
            pass
                
    def fitInView(self):
        rect = QtCore.QRectF(self.pixMap.rect())
        if not rect.isNull():
            unity = self.graphicsView.transform().mapRect(
                QtCore.QRectF(0, 0, 1, 1))
            self.graphicsView.scale(1.0 / unity.width(), 1.0 / unity.height())
            view_rect = self.graphicsView.viewport().rect()
            scene_rect = self.graphicsView.transform().mapRect(rect)
            factor = min(view_rect.width() / scene_rect.width(), 
                view_rect.height() / scene_rect.height())
            self.graphicsView.scale(factor, factor)
            self.graphicsView.centerOn(rect.center())
            self.zoom = 0
        
    def showImage(self, source_file):
        source_file = str(source_file)
        target_file = source_file.split('.')[0] + '.jpg'
    
        with wand.image.Image(filename=source_file) as img:
            img.format = 'jpeg'
            img.save(filename=target_file)
            
        img = Image.open(target_file, 'r')
        width, height = img.size
    
        self.scene.clear()
        self.zoom = 0
        
        self.imgQ = ImageQt.ImageQt(img)
        self.pixMap = QtGui.QPixmap.fromImage(self.imgQ)
        if self.pixMap and not self.pixMap.isNull():
            self.graphicsView.setDragMode(
                QtGui.QGraphicsView.ScrollHandDrag)
            self.scene.addPixmap(self.pixMap)
            self.fitInView()
        else:
            self.graphicsView.setDragMode(QtGui.QGraphicsView.NoDrag)
            self.scene.addPixmap(QtGui.QPixmap())
                
        self.scene.update()
        
        remove(target_file)


def main():
    app = QtGui.QApplication(sys.argv)
    form = PySCUBApp()
    form.show()
    sys.exit(app.exec_())
    

if __name__ == '__main__':
    main()
    
