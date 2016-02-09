# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DispersionGui.ui'
#
# Created: Mon Feb  8 09:54:39 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_DispersionCompensation(object):
    def setupUi(self, DispersionCompensation):
        DispersionCompensation.setObjectName("DispersionCompensation")
        DispersionCompensation.resize(800, 600)
        self.centralwidget = QtGui.QWidget(DispersionCompensation)
        self.centralwidget.setObjectName("centralwidget")
        self.bscan = ImageView(self.centralwidget)
        self.bscan.setGeometry(QtCore.QRect(0, 10, 611, 561))
        self.bscan.setObjectName("bscan")
        self.coef_3 = QtGui.QLineEdit(self.centralwidget)
        self.coef_3.setGeometry(QtCore.QRect(680, 10, 113, 20))
        self.coef_3.setObjectName("coef_3")
        self.coef_2 = QtGui.QLineEdit(self.centralwidget)
        self.coef_2.setGeometry(QtCore.QRect(680, 40, 113, 20))
        self.coef_2.setObjectName("coef_2")
        DispersionCompensation.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(DispersionCompensation)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 19))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        DispersionCompensation.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(DispersionCompensation)
        self.statusbar.setObjectName("statusbar")
        DispersionCompensation.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(DispersionCompensation)
        self.actionOpen.setObjectName("actionOpen")
        self.actionQuit = QtGui.QAction(DispersionCompensation)
        self.actionQuit.setObjectName("actionQuit")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionQuit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(DispersionCompensation)
        QtCore.QMetaObject.connectSlotsByName(DispersionCompensation)

    def retranslateUi(self, DispersionCompensation):
        DispersionCompensation.setWindowTitle(QtGui.QApplication.translate("DispersionCompensation", "Dispersion Compensation", None, QtGui.QApplication.UnicodeUTF8))
        self.coef_3.setText(QtGui.QApplication.translate("DispersionCompensation", "0.0e-18", None, QtGui.QApplication.UnicodeUTF8))
        self.coef_2.setText(QtGui.QApplication.translate("DispersionCompensation", "0.0e-12", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("DispersionCompensation", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("DispersionCompensation", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionQuit.setText(QtGui.QApplication.translate("DispersionCompensation", "Quit", None, QtGui.QApplication.UnicodeUTF8))

from pyqtgraph import ImageView
