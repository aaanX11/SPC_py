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

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def readConfigPechs(configFN):
    d={}
    with open(configFN) as configf:
        for line in configf:
            line = line.split('#')[0]

            if len(line)<3:
                continue
            k,v = line.strip().split('=')
            if "SPECTRS_PATH" in k:
                d["spcdir"] = v
            elif "SOURCE_FILE" in k:
                d["sourcefn"] = v
            elif "INITIALS_PATH" in k:
                d["inidir"] = v
            elif "RESULTS_PATH" in k:
                d["resdir"] = v
            elif "DETECTORS_FILE" in k:
                d["detfi"] = v
    return d
def readSourcePechs(sourceFN):
    d={}
    with open(sourceFN) as sourcef:
        for line in sourcef:
            line = line.split('#')[0]
            if len(line) < 3:
                continue
            k, v = line.strip().split('=')
            if k == "SPECTR_FILE":
               d["spcfn"] = v
    return d
def readSpcPechs(spcFN):
    d={}
    with open(spcFN) as spcf:
        while True:
            line = spcf.readline().split('#')[0]
            if len(line) < 3:
                continue
            else:
                d["type"] = line.strip().split('=')[1]
                break
        
        spcf.readline()

        d["en"] = []
        d["fr"] = []
        if "CONTINUOUS" == d["type"]:
            d["en"].append(float(spcf.readline().strip()))
        #if "DISCRETE" == d["type"]:
        lines = spcf.readlines()
        for line in lines:
            en, fr = map(lambda x: float(x.strip()), line.split())
            d["en"].append(en)
            d["fr"].append(fr)
    return d

def read_dets(path):
    dets = []
    with open(path, 'r') as detf:
        d = {}
        for line in detf:
            line = line.split('#')[0].strip()

            if len(line)<3:
                dets.append(d)
                continue
            if '[FAMILY]' in line:
                d = {}
                continue
            k,v = line.strip().split('=')
            if "FAMILY_NAME" in k:
                d["name"] = v
            elif "PARTICLE_TYPE" in k:
                d["particle"] = v
            elif "DETECTOR_TYPE" in k:
                d["type"] = v
            elif "DETECTOR_SHAPE_TYPE" in k:
                d["shape"] = v
            elif "FLUX_MEASUREMENT" in k:
                d["measure"] = v
            elif "FLUX_DIRECTIONS" in k:
                d["directn"] = v
            elif "LIST_FILE" in k:
                d["list"] = v
            elif "TEMPLATE_FILE" in k:
                d["templ"] = v
    return dets                

def readPECHS(path):

    config = readConfigPechs(os.path.join(path, "configuration"))
    
    source = readSourcePechs(os.path.join(path,  config["inidir"], config["sourcefn"]))
    
    spc = {}#readSpcPechs(os.path.join(path, config["spcdir"], source["spcfn"]))

    dets = read_dets(os.path.join(path, config["inidir"], config["detfi"] ))

    return config, source, spc, dets
    
