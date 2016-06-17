#!/usr/bin/env python


# PySCUBA/src/PySCUBA/__main__.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from os import getcwd, path
import sys

try:
    import igraph
except ImportError, e:
    pass

import numpy as np
from PyQt4 import QtCore, QtGui
from sklearn.preprocessing import StandardScaler

from .Gap_stats import gap_stats
from .Preprocessing import cytometry_preprocess, PCR_preprocess, RNASeq_preprocess
import .PySCUBA_design
import .SCUBA_core as SCUBA


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

    def __init__(self, data_type, data_path, cluster_mode, log_mode,
                 pseudotime_mode, pcv_method, anchor_gene,
                 exclude_marker_names):
        super(WorkerThread, self).__init__()
        
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
            PCA_components, data = SCUBA.PCA_analysis(data, self.cluster_mode,
                cell_stages if (self.cluster_mode == 'pca2') else None)
    
        centroid_coords, cluster_indices, \
        parent_clusters = SCUBA.initialize_tree(data, cell_stages)
        centroid_coords, cluster_indices, \
        parent_clusters = SCUBA.refine_tree(data, centroid_coords,
            cluster_indices, parent_clusters, cell_stages, output_directory)
    
        plot_tree(cluster_indices, parent_clusters, output_directory)
        
        if self.cluster_mode in {'pca', 'pca2'}:
            weights = PCA_components
        else:
            weights = None
    
        bifurcation_info, bifurcation_axes, \
        bifurcation_projections = SCUBA.bifurcation_direction(data, cell_IDs,
            markers, parent_clusters, centroid_coords, output_directory,
            weights)
                
        if bifurcation_info:
            data_per_split, parameters_per_split = SCUBA.bifurcation_analysis(
                cluster_indices, bifurcation_info, bifurcation_axes,
                bifurcation_projections, output_directory)
        
        return


class PySCUBApp(QtGui.QMainWindow, PySCUBA_design.Ui_MainWindow):

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)

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
        
        self.okButton.clicked.connect(self.buttonClicked)
        self.okButton.clicked.connect(self.OK)
        self.cancelButton.clicked.connect(self.buttonClicked)
        
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
        self.dataFileDialog = QtGui.QFileDialog(self)
        self.data_path = str(self.dataFileDialog.getOpenFileName())
        self.statusbar.showMessage("{0} ready to be "
            "analyzed".format(path.basename(self.data_path)))
    
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
        
        self.get_thread = WorkerThread(self.dataTypeComboBox.currentText(),
            self.data_path, self.clusterModeComboBox.currentText(),
            self.log_mode, self.pseudotime_mode, self.pcv_method,
            self.anchor_gene, self.exclude_marker_names)
        self.connect(self.get_thread, QtCore.SIGNAL("update(QString)"),
            self.get_thread.run)
        self.connect(self.get_thread, QtCore.SIGNAL("finished()"), self.done)
        self.get_thread.start()
        self.cancelButton.setEnabled(True)
        self.okButton.setEnabled(False)
    
    def done(self):
        if self.button_clicked == 'Cancel':
            self.cancel()
        else:
            self.cancelButton.setEnabled(False)
            self.okButton.setEnabled(True)
            self.statusbar.showMessage("PySCUBA has completed the "
                "analysis of your data.")
            QtGui.QMessageBox.information(self, "Status Message", 
                "Mission completed!")
            
    def cancel(self):
        self.get_thread.terminate()
        self.statusbar.showMessage("PySCUBA was interrupted!")
        QtGui.QMessageBox.information(self, "Status Message", 
            "PySCUBA was interrupted!")
    
    def selectDisplay(self):
        pass


def main():
    app = QtGui.QApplication(sys.argv)
    form = PySCUBApp()
    form.show()
    sys.exit(app.exec_())
    

if __name__ == '__main__':
    main()
    
