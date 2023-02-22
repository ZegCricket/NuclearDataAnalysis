import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from pybaselines.smooth import snip
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys
from sasnip import sasnip

def calculate_snip_background(data, mhw, shw, decrease = True, lls = True, fo = 2):
    if lls == True:
        transformed_data =  np.log(np.log(np.sqrt(data + 1) + 1) + 1)
        background = snip(transformed_data, max_half_window = mhw, decreasing = decrease, smooth_half_window = shw, filter_order = fo)[0]
        background = -1 + (np.exp(np.exp(background) - 1) - 1) ** 2
    else:
        background = snip(data, max_half_window = mhw, decreasing = decrease, smooth_half_window = shw, filter_order = fo)[0]
    data_no_bkg = data - background
    return data_no_bkg, background

def onlyFileName(fileName):
    stringSize = len(fileName)
    for i in range(stringSize):
        if fileName[stringSize - (i + 1)] == '/':
            dirName = fileName[:stringSize - (i + 1)]
            fileName = fileName[stringSize - i:]
            break
    return fileName, dirName

def twoSpinBoxes(labels, layout):
    label1 = QLabel(labels[0])
    label1.setFixedHeight(15)
    label1.setAlignment(Qt.AlignRight)
    layout.addWidget(label1)
    spinBox1 = QSpinBox()
    spinBox1.setFixedWidth(75)
    spinBox1.setAlignment(Qt.AlignCenter)
    layout.addWidget(spinBox1)
    label2 = QLabel(labels[1])
    label2.setFixedHeight(15)
    label2.setAlignment(Qt.AlignRight)
    layout.addWidget(label2)
    spinBox2 = QSpinBox()
    spinBox2.setFixedWidth(75)
    spinBox2.setAlignment(Qt.AlignCenter)
    layout.addWidget(spinBox2)
    return spinBox1, spinBox2

class window(QMainWindow):
    def __init__(self, parent = None):
        super(window, self).__init__(parent)

        #Main window setup
        win = QWidget(self)
        self.setGeometry(100, 100, 1280, 720)
        self.span = []
        self.facecolors = plt.colormaps['winter'](np.linspace(0, 1, 4))
        self.dirName = ""
        #self.setStyleSheet("background-color: white;")

        #Menu bar setup
        menu = self.menuBar()
        #File menu
        file = menu.addMenu("File")
        open = file.addAction("Open")
        file.addAction("Save")

        #Central widget
        self.setCentralWidget(win)
        hbox = QHBoxLayout()

        #Configuration of the left part of the central widget (a.k.a. the plot)
        plotLayout = QVBoxLayout()
        hbox.addLayout(plotLayout)
        self.static_canvas = FigureCanvas(plt.Figure())
        self.navigation = NavigationToolbar(self.static_canvas)
        plotLayout.addWidget(self.navigation)
        plotLayout.addWidget(self.static_canvas)
        self.ax = self.static_canvas.figure.subplots()
        self.ax.grid(visible = True, which = "major", axis = "both")
        self.ax.set_xlabel("Channel")
        self.ax.set_ylabel("Counts")
        self.plot, = self.ax.plot(0, color = self.facecolors[0])
        self.sniplot, = self.ax.plot(0, "--", label = "SNIP", color = self.facecolors[2])
        self.nobkgplot, = self.ax.plot(0, label = "W/o Background", color = self.facecolors[3])

        #Configuration of the right part of the central widget
        toolsLayout = QVBoxLayout()
        hbox.addLayout(toolsLayout)
        #Configuration of the combo box
        self.cb = QComboBox(win)
        self.cb.setFixedWidth(300)
        self.cb.setEditable(True)
        self.cb.lineEdit().setAlignment(Qt.AlignCenter)
        self.cb.lineEdit().setReadOnly(True)
        toolsLayout.addWidget(self.cb)
        #Configuration of the "counts part"
        countsLayout = QVBoxLayout()
        toolsLayout.addLayout(countsLayout)
        countsLine = QHBoxLayout()
        countsLayout.addLayout(countsLine)
        countsButton = QPushButton("Count")
        countsButton.setFixedWidth(100)
        countsLine.addWidget(countsButton)
        self.countsString = QLineEdit("0")
        self.countsString.setFixedWidth(100)
        self.countsString.setReadOnly(True)
        self.countsString.setAlignment(Qt.AlignCenter)
        countsLine.addWidget(self.countsString)
        #Configuration of the count interval
        intervalLabel = QLabel("Count Interval")
        intervalLabel.setFixedHeight(20)
        intervalLabel.setAlignment(Qt.AlignCenter)
        countsLayout.addWidget(intervalLabel)
        intervalLine = QHBoxLayout()
        countsLayout.addLayout(intervalLine)
        minMaxLabels = np.array(["Minimum", "Maximum"])
        self.countInterval = twoSpinBoxes(minMaxLabels, intervalLine)
        #Radio buttons
        radioLayout = QHBoxLayout()
        countsLayout.addLayout(radioLayout)
        self.originalButton = QRadioButton("Original")
        self.snipButton = QRadioButton("SNIP")
        self.sasnipButton = QRadioButton("SASNIP")
        radioLayout.addWidget(self.originalButton)
        radioLayout.addWidget(self.snipButton)
        radioLayout.addWidget(self.sasnipButton)
        self.sasnipButton.setChecked(True)

        #Configuration of the SNIP options
        snipLayout = QVBoxLayout()
        toolsLayout.addLayout(snipLayout)
        snipLabel = QLabel("SNIP Options")
        snipLabel.setFixedHeight(20)
        snipLabel.setAlignment(Qt.AlignCenter)
        snipLayout.addWidget(snipLabel)
        checkBoxesLabels = np.array(["SNIP Graphic", "SNIP Count", "Decreased", "LLS Operator"])
        self.checkBoxes = np.zeros(checkBoxesLabels.size, dtype = QCheckBox)
        for i in range(checkBoxesLabels.size):
            self.checkBoxes[i] = QCheckBox(win)
            self.checkBoxes[i].setText(checkBoxesLabels[i])
            if i != 1:
                self.checkBoxes[i].setChecked(True)
            snipLayout.addWidget(self.checkBoxes[i])
        #Configuration of the SNIP interval
        snipIntervalLayout = QHBoxLayout()
        snipLayout.addLayout(snipIntervalLayout)
        twoSpinBoxes(minMaxLabels, snipIntervalLayout)
        #Configuration of the SNIP parameters
        snipParameterLine = QHBoxLayout()
        snipLayout.addLayout(snipParameterLine)
        snipParameterLabels = np.array(["MHWindow", "Smooth"])
        self.mhWindow, self.smooth = twoSpinBoxes(snipParameterLabels, snipParameterLine)
        self.mhWindow.setMinimum(1)
        self.mhWindow.setValue(10)
        self.smooth.setValue(3)

        sasnipLayout = QVBoxLayout()
        toolsLayout.addLayout(sasnipLayout)
        sasnipLabel = QLabel("SASNIP Options")
        sasnipLabel.setFixedHeight(20)
        sasnipLabel.setAlignment(Qt.AlignCenter)
        sasnipLayout.addWidget(sasnipLabel)
        sacheckBoxesLabels = np.array(["SASNIP Graphic", "Decreased"])
        self.sacheckBoxes = np.zeros(sacheckBoxesLabels.size, dtype = QCheckBox)
        for i in range(sacheckBoxesLabels.size):
            self.sacheckBoxes[i] = QCheckBox(win)
            self.sacheckBoxes[i].setText(sacheckBoxesLabels[i])
            self.sacheckBoxes[i].setChecked(True)
            sasnipLayout.addWidget(self.sacheckBoxes[i])
        #Configuration of the SASNIP parameters
        sasnipParameterLine = QHBoxLayout()
        sasnipLayout.addLayout(sasnipParameterLine)
        sasnipParameterLabels = np.array(["Max", "Threshold"])
        self.max, self.threshold = twoSpinBoxes(sasnipParameterLabels, sasnipParameterLine)
        self.max.setMinimum(1)
        self.max.setMaximum(200)
        self.max.setValue(100)
        self.threshold.setValue(1)

        #Triggers
        open.triggered.connect(self.openFile)
        countsButton.clicked.connect(self.count)
        self.countInterval[0].valueChanged.connect(self.updateROI)
        self.countInterval[1].valueChanged.connect(self.updateROI)
        self.mhWindow.valueChanged.connect(self.updateSNIP)
        self.smooth.valueChanged.connect(self.updateSNIP)
        self.checkBoxes[2].stateChanged.connect(self.updateSNIP)
        self.checkBoxes[3].stateChanged.connect(self.updateSNIP)
        self.originalButton.toggled.connect(self.updateSNIP)
        self.snipButton.toggled.connect(self.updateSNIP)
        self.sasnipButton.toggled.connect(self.updateSNIP)
        self.max.valueChanged.connect(self.updateSNIP)
        self.threshold.valueChanged.connect(self.updateSNIP)
        self.sacheckBoxes[1].stateChanged.connect(self.updateSNIP)

        #Config
        win.setLayout(hbox)
        self.setWindowTitle("Nuclear Data Analysis")
        self.showMaximized()
    
    def openFile(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileNames(self, "Open File", self.dirName, "All Files (*);;Text Files (*.txt);;Data Files (*.dat)", options = options)
        if fileName[0] != "":
            try:
                self.data = np.genfromtxt(fileName[0], comments = '$')
            except ValueError:
                self.data = np.genfromtxt(fileName[0], comments = '$', skip_footer = 1)
            self.data = np.reshape(self.data, -1)
            fileName, self.dirName = onlyFileName(fileName[0])
            if self.sasnipButton.isChecked():
                background = sasnip(self.data, peakMaximum = self.max.value(), derivativeThreshold = self.threshold.value(), decrease = self.sacheckBoxes[1].isChecked())
                self.data_no_bkg = self.data - background
            else:
                self.data_no_bkg, background = calculate_snip_background(self.data, self.mhWindow.value(), self.smooth.value(), self.checkBoxes[2].isChecked(), self.checkBoxes[3].isChecked())
            x = range(0, self.data.size)
            self.plot.set_data(x, self.data)
            self.plot.set_label(fileName)
            self.sniplot.set_data(x, background)
            self.nobkgplot.set_data(x, self.data_no_bkg)
            self.updateROI()
            self.navigation.update()
            self.ax.set_xlim(0, self.data.size)
            self.ax.set_ylim(0, self.data.max() * 1.1)
            self.ax.legend()
            self.static_canvas.draw_idle()
            self.countInterval[0].setMaximum(self.data.size - 1)
            self.countInterval[1].setMaximum(self.data.size - 1)
            self.cb.clear()
            self.cb.addItem(fileName)

    def count(self):
        if self.cb.count() != 0:
            area = 0
            if (self.snipButton.isChecked() or self.sasnipButton.isChecked()):
                for i in range(self.countInterval[0].value(), self.countInterval[1].value(), 1):
                    area += self.data_no_bkg[i]
            else:
                for i in range(self.countInterval[0].value(), self.countInterval[1].value(), 1):
                    area += self.data[i]
            self.countsString.setText(str(int(area)))
    
    def updateSNIP(self):
        if self.cb.count() != 0:
            x = range(0, self.data.size)
            if self.sasnipButton.isChecked():
                background = sasnip(self.data, peakMaximum = self.max.value(), derivativeThreshold = self.threshold.value(), decrease = self.sacheckBoxes[1].isChecked())
                self.data_no_bkg = self.data - background
            else:
                self.data_no_bkg, background = calculate_snip_background(self.data, self.mhWindow.value(), self.smooth.value(), self.checkBoxes[2].isChecked(), self.checkBoxes[3].isChecked())
            self.sniplot.set_data(x, background)
            self.nobkgplot.set_data(x, self.data_no_bkg)
            self.static_canvas.draw_idle()

    def updateROI(self):
        if len(self.span):
            for i in self.span:
                try:
                    i.remove()
                except ValueError:
                    continue
        self.span.append(self.ax.axvspan(self.countInterval[0].value(), self.countInterval[1].value(), alpha = 0.5, color = self.facecolors[1]))
        self.static_canvas.draw_idle()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = window()
    win.show()
    sys.exit(app.exec_())