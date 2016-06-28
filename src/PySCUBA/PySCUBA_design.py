#!/usr/bin/env python
# -*- coding: utf-8 -*-


# PySCUBA/src/PySCUBA/PySCUBA_design.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from PyQt4 import QtCore, QtGui


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, 
            disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)


class Ui_MainWindow(object):

    def setupUi(self, MainWindow):
        self.setupMainWindow(MainWindow)
        
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        
        self.setupMenuBar(MainWindow)
        self.setupStatusBar(MainWindow)
        
        self.mainVerticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.mainVerticalLayout.setObjectName(_fromUtf8("mainVerticalLayout"))
        
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setHorizontalSpacing(3)
        self.gridLayout.setVerticalSpacing(2)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        
        self.hline1 = QtGui.QFrame(self.centralwidget)
        self.hline1.setFrameShadow(QtGui.QFrame.Raised)
        self.hline1.setLineWidth(7)
        self.hline1.setFrameShape(QtGui.QFrame.HLine)
        self.hline1.setObjectName(_fromUtf8("hline1"))
        self.hline1.raise_()
        
        self.cancelOkLayout = QtGui.QHBoxLayout()
        self.cancelOkLayout.setObjectName(_fromUtf8("cancelOkLayout"))
        
        self.hline2 = QtGui.QFrame(self.centralwidget)
        self.hline2.setFrameShadow(QtGui.QFrame.Raised)
        self.hline2.setLineWidth(7)
        self.hline2.setFrameShape(QtGui.QFrame.HLine)
        self.hline2.setObjectName(_fromUtf8("hline2"))
        self.hline2.raise_()
        
        self.graphicsVerticalLayout = QtGui.QVBoxLayout()
        self.graphicsVerticalLayout.setObjectName(
            _fromUtf8("graphicsVerticalLayout"))
        
        self.mainVerticalLayout.addLayout(self.gridLayout)
        self.mainVerticalLayout.addWidget(self.hline1)
        self.mainVerticalLayout.addLayout(self.cancelOkLayout)
        self.mainVerticalLayout.addWidget(self.hline2)
        self.mainVerticalLayout.addLayout(self.graphicsVerticalLayout)
        
        self.adornGridLayout(MainWindow)
        self.adornCancelOkLayout(MainWindow)
        self.adornGraphicsVerticalLayout(MainWindow)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
    def setupMainWindow(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.setWindowModality(QtCore.Qt.NonModal)
        MainWindow.setGeometry(150, 100, 564, 635)
        MainWindow.setMouseTracking(False)
        MainWindow.setFocusPolicy(QtCore.Qt.ClickFocus)
        MainWindow.setAutoFillBackground(False)
        MainWindow.setTabShape(QtGui.QTabWidget.Rounded)
        
    def setupMenuBar(self, MainWindow):
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionExit.setShortcut('Ctrl+Q')
        self.actionExit.setStatusTip('Exit application')
        self.actionExit.triggered.connect(QtGui.qApp.quit)
        
        self.actionCredits = QtGui.QAction(MainWindow)
        self.actionCredits.setObjectName(_fromUtf8("actionCredits"))
        self.actionCredits.setShortcut('Ctrl+I')
        self.actionCredits.setStatusTip('Display credits')
        self.actionCredits.triggered.connect(self.showCredits)
        
        self.actionHelp = QtGui.QAction(MainWindow)
        self.actionHelp.setObjectName(_fromUtf8("actionHelp"))
        self.actionHelp.setShortcut('Ctrl+H')
        self.actionHelp.setStatusTip('Help and documentation')
        self.actionHelp.triggered.connect(self.showDocumentation)
        
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.fileMenu = self.menubar.addMenu('&File')
        self.fileMenu.addAction(self.actionExit)
        self.fileMenu = self.menubar.addMenu('&Credits')
        self.fileMenu.addAction(self.actionCredits)
        self.fileMenu = self.menubar.addMenu('&Help')
        self.fileMenu.addAction(self.actionHelp)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 542, 23))
        self.menubar.setDefaultUp(False)
        self.menubar.setNativeMenuBar(False)
        self.menubar.setObjectName(_fromUtf8("menubar"))
        
        MainWindow.setMenuBar(self.menubar)
    
    def setupStatusBar(self, MainWindow):
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.showMessage(
            "Ready - Please select a dataset to process")
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        
        MainWindow.setStatusBar(self.statusbar)
      
    def showCredits(self, MainWindow):
        QtGui.QMessageBox.information(self, "Credits", 
            "Author: Gregory Giecold\nAffiliation: Harvard "
            "University & DFCI\nContact: ggiecold@jimmy.harvard.edu\n"
            "GitHub: https://github.com/GGiecold\n")
            
    def showDocumentation(self, MainWindow):
        url = QtCore.QUrl("https://github.com/GGiecold/PySCUBA")
        if not QtGui.QDesktopServices.openUrl(url):
            QtGui.QMessageBox.warning(self, 'Help & Documentation', 
                'Could not open url to online documentation!')
            
    def adornGridLayout(self, MainWindow):
        self.datasetHorizontalLayout = QtGui.QHBoxLayout()
        self.datasetHorizontalLayout.setObjectName(
            _fromUtf8("datasetHorizontalLayout"))
        self.gridLayout.addLayout(self.datasetHorizontalLayout, 0, 0, 1, 1)
        
        self.datasetHorizontalLayout.addStretch(1)
        self.selectDatasetButton = QtGui.QPushButton(self.centralwidget)
        self.selectDatasetButton.setToolTip("Click this button to browse "
            "through\nyour computer's directories and select\na dataset "
            "to subject to a clustering\nand bifurcation analysis.")
        self.selectDatasetButton.setObjectName(
            _fromUtf8("selectDatasetButton"))
        self.datasetHorizontalLayout.addWidget(self.selectDatasetButton)
        self.datasetHorizontalLayout.addStretch(1)
        
        self.withinGridVerticalLayout_1 = QtGui.QVBoxLayout()
        self.withinGridVerticalLayout_1.setObjectName(
            _fromUtf8("withinGridVerticalLayout_1"))
        self.gridLayout.addLayout(self.withinGridVerticalLayout_1, 1, 0, 1, 1)
            
        self.dataTypeLabel = QtGui.QLabel(self.centralwidget)
        self.dataTypeLabel.setFrameShadow(QtGui.QFrame.Raised)
        self.dataTypeLabel.setLineWidth(5)
        self.dataTypeLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.dataTypeLabel.setObjectName(_fromUtf8("dataTypeLabel"))
        self.withinGridVerticalLayout_1.addWidget(self.dataTypeLabel)
        
        self.dataTypeComboBox = QtGui.QComboBox(self.centralwidget)
        self.dataTypeComboBox.setToolTip("Does the file to process qualify "
            "as qPCR, RNAseq\nor flow or mass cytometry data?\nThe latter "
            "is expected to be in *.fcs format,\nwhereas the first two "
            "types should be delivered\nas rows of tab-separated entries.")
        self.dataTypeComboBox.setEditable(True)
        self.dataTypeComboBox.setObjectName(_fromUtf8("dataTypeComboBox"))
        self.dataTypeComboBox.addItem(_fromUtf8(""))
        self.dataTypeComboBox.addItem(_fromUtf8(""))
        self.dataTypeComboBox.addItem(_fromUtf8(""))
        self.withinGridVerticalLayout_1.addWidget(self.dataTypeComboBox)
        
        self.withinGridVerticalLayout_2 = QtGui.QVBoxLayout()
        self.withinGridVerticalLayout_2.setObjectName(
            _fromUtf8("withinGridVerticalLayout_2"))
        self.gridLayout.addLayout(self.withinGridVerticalLayout_2, 2, 0, 1, 1)
        
        self.clusterModeLabel = QtGui.QLabel(self.centralwidget)
        self.clusterModeLabel.setFrameShadow(QtGui.QFrame.Raised)
        self.clusterModeLabel.setLineWidth(5)
        self.clusterModeLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.clusterModeLabel.setObjectName(_fromUtf8("clusterModeLabel"))
        self.withinGridVerticalLayout_2.addWidget(self.clusterModeLabel)
        
        self.clusterModeComboBox = QtGui.QComboBox(self.centralwidget)
        self.clusterModeComboBox.setToolTip("For each timestamp binned to "
            "a particular stage, PySCUBA proceeds\nto several rounds of "
            "aggregating the corresponding samples\ninto an optimal number "
            "of clusters.\nBy selecting 'PCA' from this menu, the clustering "
            "will be based on\na reduction of the original dataset to its "
            "first few principal components.\nA choice of 'PCA2' specifies "
            "that the principal components analysis\nwill be based on "
            "the samples that are part of the final stage\nof the "
            "temporal ordering.")
        self.clusterModeComboBox.setEditable(True)
        self.clusterModeComboBox.setObjectName(
            _fromUtf8("clusterModeComboBox"))
        self.clusterModeComboBox.addItem(_fromUtf8(""))
        self.clusterModeComboBox.addItem(_fromUtf8(""))
        self.clusterModeComboBox.addItem(_fromUtf8(""))
        self.withinGridVerticalLayout_2.addWidget(self.clusterModeComboBox)
        
        self.logCheckBox = QtGui.QCheckBox(self.centralwidget)
        self.logCheckBox.setChecked(True)
        self.logCheckBox.setObjectName(_fromUtf8("logCheckBox"))
        self.gridLayout.addWidget(self.logCheckBox, 0, 1, 1, 1)
        
        self.pseudotimeCheckBox = QtGui.QCheckBox(self.centralwidget)
        self.pseudotimeCheckBox.setToolTip("If your data is not endowed with "
            "temporal information of any kind, please\ndo check this box. "
            "PySCUBA will thereby run a principal curve analysis\nto infer a "
            "temporal ordering for each sample of your dataset.")
        self.pseudotimeCheckBox.setChecked(True)
        self.pseudotimeCheckBox.setObjectName(_fromUtf8("pseudotimeCheckBox"))
        self.gridLayout.addWidget(self.pseudotimeCheckBox, 1, 1, 1, 1)

    def adornCancelOkLayout(self, MainWindow):
        self.cancelOkLayout.addStretch(1)
        self.cancelButton = QtGui.QPushButton(self.centralwidget)
        self.cancelButton.setObjectName(_fromUtf8("cancelButton"))
        self.cancelOkLayout.addWidget(self.cancelButton)
        
        self.okButton = QtGui.QPushButton(self.centralwidget)
        self.okButton.setToolTip("Click this button to browse "
            "through\nyour computer's directories and select\na dataset "
            "to subject to a clustering\nand bifurcation analysis.")
        self.okButton.setObjectName(_fromUtf8("okButton"))
        self.cancelOkLayout.addWidget(self.okButton)
        self.cancelOkLayout.addStretch(1)
        
    def adornGraphicsVerticalLayout(self, MainWindow):
        self.scene = QtGui.QGraphicsScene(self.centralwidget)

        self.graphicsView = QtGui.QGraphicsView(self.scene)
        self.graphicsView.setFrameShadow(QtGui.QFrame.Raised)
        self.graphicsView.setLineWidth(3)
        self.graphicsView.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOn)
        self.graphicsView.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOn)
        self.graphicsView.setTransformationAnchor(
            QtGui.QGraphicsView.AnchorUnderMouse)
        self.graphicsView.setResizeAnchor(
            QtGui.QGraphicsView.AnchorUnderMouse)
        self.graphicsView.setBackgroundBrush(
            QtGui.QBrush(QtGui.QColor(245,245,245)))
        self.graphicsView.setFrameShape(QtGui.QFrame.NoFrame)
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.graphicsVerticalLayout.addWidget(self.graphicsView)
        
        self.displayHorizontalLayout = QtGui.QHBoxLayout()
        self.displayHorizontalLayout.setObjectName(
            _fromUtf8("displayHorizontalLayout"))
        self.graphicsVerticalLayout.addLayout(self.displayHorizontalLayout)
        
        self.displayHorizontalLayout.addStretch(1)
        self.displayFileButton = QtGui.QPushButton(self.centralwidget)
        self.displayFileButton.setToolTip("Various files and figures will "
            "show up in this box as they are\nbeing produced by the PySCUBA "
            "analysis of your data.\nClick on any of those and it will be "
            "displayed in an adjacent\ngraphics box.")
        self.displayFileButton.setObjectName(_fromUtf8("displayFileButton"))
        self.displayHorizontalLayout.addWidget(self.displayFileButton)
        self.displayHorizontalLayout.addStretch(1)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", 
            "PySCUBA - GC Yuan Lab", None))
            
        self.selectDatasetButton.setText(_translate("MainWindow", 
            "1. Select dataset to analyze", None))
        
        self.dataTypeLabel.setText(_translate("MainWindow", 
            "2. Specify type of data:", None))
        self.dataTypeComboBox.setItemText(0, _translate("MainWindow",
            "RNASeq", None))
        self.dataTypeComboBox.setItemText(1, _translate("MainWindow", 
            "PCR", None))
        self.dataTypeComboBox.setItemText(2, _translate("MainWindow",
            "cytometry", None))
        
        self.clusterModeLabel.setText(_translate("MainWindow", 
            "3. Choose cluster mode:", None))
        self.clusterModeComboBox.setItemText(0, _translate("MainWindow",
            "None", None))
        self.clusterModeComboBox.setItemText(1, _translate("MainWindow",
            "PCA", None))
        self.clusterModeComboBox.setItemText(2, _translate("MainWindow",
            "PCA2", None))
        
        self.logCheckBox.setText(_translate("MainWindow", 
            "4. Apply a log-transform?", None))
        
        self.pseudotimeCheckBox.setText(_translate("MainWindow", 
            "5. Infer temporal ordering?", None))
            
        self.cancelButton.setText(_translate("MainWindow", "Cancel", None))
        self.okButton.setText(_translate("MainWindow", "Ok", None))
        
        self.displayFileButton.setText(_translate("MainWindow", "Select file to display", None))
        
        self.actionExit.setText(_translate("MainWindow", "Exit", None))


if __name__ == "__main__":
    import sys
    
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())


