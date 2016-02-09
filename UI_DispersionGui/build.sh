#! /bin/bash

rm UI_DispersionGui.py
pyside-uic DispersionGui.ui > UI_DispersionGui.py
cp UI_DispersionGui.py ..
