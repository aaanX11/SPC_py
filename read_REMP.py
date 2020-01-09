import numpy as np

import matplotlib.pyplot as plt

from matplotlib.widgets import Button
from matplotlib.widgets import TextBox

from scipy import integrate

from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

import sys

import os
from Tkinter import *
import tkMessageBox

def is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_start(startFN):
    d={}
    with open(startFN) as startf:
        d["prjf"] = startf.readline().strip()
        # d["npx"] = int(startf.readline().strip())
    return d


def read_prj(prjFN):
    d={}
    with open(prjFN) as prjf:
        
        line = prjf.readline().strip()
        while line:
            line = prjf.readline().strip()
            # print line.strip("><")
            if "Dxf name" in line:
                d["dxff"] = prjf.readline().strip()
            elif "Grd name" in line:
                d["grdf"] = prjf.readline().strip()
            elif "Cell name" in line:
                d["celf"] = prjf.readline().strip()
            elif "Particles-Layers name" in line:
                d["plf"] = prjf.readline().strip()
            elif "Layers name" in line:
                d["layf"] = prjf.readline().strip()
            elif "Output" in line:
                n_outp = int(prjf.readline().strip())
                if n_outp > 0:
                    d["vidf"] = []
                for i in range(n_outp):
                    d["vidf"].append(prjf.readline().strip())
                    prjf.readline()
    return d


def read_grid(grdFN):
    d = {}
    with open(grdFN) as grdf:
        grdf.readline()  # <FULL>
        grdf.readline()  # 1
        grdf.readline()  # <EQUAL (0-not equal 1-equal)>
        grdf.readline()  # 0
        grdf.readline()  # X
        d["nx"] = int(grdf.readline())
        d["x"] = map(float, grdf.readline().strip().split())

        grdf.readline()
        d["ny"] = int(grdf.readline())
        d["y"] = map(float, grdf.readline().strip().split())
        
        grdf.readline()
        d["nz"] = int(grdf.readline())
        d["z"] = map(float, grdf.readline().strip().split())

        grdf.readline()
        d["nt"] = int(grdf.readline())
        d["t"] = map(float, grdf.readline().strip().split())
    return d


def process_grid(xi):
    
    d = {}
    
    xi05 = [0.5*(x1+x2) for x1,x2 in zip(xi[1:], xi[:-1])]

    dxi05 = [x2-x1 for x2,x1 in zip(xi[1:], xi[:-1])]
    dxi05.insert(0,dxi05[0])
    dxi05.append(dxi05[-1])
    
    xi05.insert(0, xi[0]-0.5*dxi05[0])
    xi05.append(xi[-1]+0.5*dxi05[-1])
    d['i'] = xi
    d['i05'] = xi05
    d['di05'] = dxi05
    
    return d


def spaceOf(celFN, shape):
    cells = np.loadtxt(celFN, np.dtype(int))
    size = (shape[0]-2)*(shape[1]-2)*(shape[2]-2)
    space1 = np.array([-1]*size)
    space = np.full(shape, -1)
    
    assert np.sum(cells[1::2]) == size
   
    filled = 0
    for i in range(cells.size/2):
        fill = cells[2*i+1]
        space1[filled:filled+fill] = cells[2*i]
        filled += fill
    space1.resize((shape[0]-2), (shape[1]-2), (shape[2]-2))
    space[1:-1, 1:-1, 1:-1] = space1
    return space


def readREMP(path):
    if os.path.exists(os.path.join(path, "START_N")):
        start = read_start(os.path.join(path, "START_N"))
    else:
        try:
            start = read_start(os.path.join(path, "START"))
        except (IOError, ValueError) as e: 
            print repr(e)
            exit()
    prj = read_prj(os.path.join(path, start["prjf"]))
    grd = read_grid(os.path.join(path, prj["grdf"]))
    prj['prjf'] = start["prjf"]
    space = spaceOf(os.path.join(path, prj["celf"]), (grd['nx']+2, grd['ny']+2, grd['nz']+2))
    return prj, grd, space


    
