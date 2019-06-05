import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize

import Tkinter as tk
import tkMessageBox
import ttk



import matplotlib
import matplotlib.pyplot as plt

from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import matplotlib.collections as clt
import mpl_toolkits.mplot3d as a3

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import read_REMP
import read_PECHS
import spc_en
import spc_tok

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_dets(path):
    dets = []
    with open(path, 'r') as detf:
        d = {}
        for line in detf:
            line = line.split('#')[0].strip()

            if len(line) < 3:
                dets.append(d)
                continue
            if '[FAMILY]' in line:
                d = {}
                continue
            k, v = line.strip().split('=')
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

def read_tmpl(path):
    with open(path, 'r') as f:
        d = {}
        for line in f:
            line = line.split('#')[0].strip()

            if len(line) < 3:
                break
            k, v = line.strip().split('=')
            if "ENERGY_GRID_TYPE" in k:
                d["type"] = v
            elif "UNIFORM_MIN_ENERGY" in k:
                d["minen"] = float(v)
            elif "UNIFORM_MAX_ENERGY" in k:
                d["maxen"] = float(v)
            elif "UNIFORM_ENERGY_RES" in k:
                d["res"] = int(v)
            elif "CUSTOM" in k:
                d["energies"] = map(float, v.split())
    if 'energies' not in d:
        interv = (d['maxen'] - d['minen']) / d['res']
        d['energies'] = [d['minen'] + 0.5 * interv + i * interv for i in range(d['res'])]
    return d
    

def read_det_info(path, config, d):
    lstf = np.loadtxt(os.path.join(path, config["inidir"], d['list']))
    totf = np.loadtxt(os.path.join(path, config["resdir"], d['name']+'.tot'))
    resf = np.genfromtxt(os.path.join(path, config["resdir"], d['name']+'.res'))
    varf = np.loadtxt(os.path.join(path, config["resdir"], d['name']+'.var'))
    return lstf, totf, resf, varf


def count_directns_b(lstf):
    directns = np.zeros((3, 2))
    maskx = (abs(lstf[:, 3]) > 1e-32)
    # = np.unique(maskx, return_counts = True)
    _, ind = np.unique(maskx, return_counts=True)
    directns[0, :] = ind
    masky = (abs(lstf[:, 4]) > 1e-32)
    _, ind = np.unique(masky, return_counts=True)
    directns[1, :] = ind
    print directns[1, :]
    maskz = (abs(lstf[:, 5]) > 1e-32)
    _, ind = np.unique(maskz, return_counts=True)
    directns[2, :] = ind
    print directns[2, :]
    return directns

def count_directns_d(resf):
    directns = np.zeros((3, 2))
    maskx = (abs(resf[:, 1]) > 1e-32)
    # = np.unique(maskx, return_counts = True)
    _, ind = np.unique(maskx, return_counts=True)
    directns[0, :] = ind
    masky = (abs(resf[:, 2]) > 1e-32)
    _, ind = np.unique(maskx, return_counts=True)
    directns[1, :] = ind
    print directns[1,:]
    maskz = (abs(resf[:, 3]) > 1e-32)
    _, ind = np.unique(maskz, return_counts=True)
    directns[2, :] = ind
    print directns[2, :]
    return directns

def process_b_fl(name, tmpl, directns, lstf, totf, resf, varf):
    energies = tmpl['energies']

    ptcl = 0
    files = {}
    for dir_, dir_cnt in enumerate(directns.flatten()):
        if dir_cnt > 0:
            files[dir_] = open(name+'_' +str(dir_)+'_'+ str(ptcl)+".spc", 'w')
            print name +str(dir_)+'_'+ str(ptcl)+".spc"

    for listi, toti, resi, vari in zip(lstf, totf, resf, varf):
        if len(list(resi)) < len(energies):
            assert('Number of entries in detector is less than expected. Check detector type')
            
        entry = '\t'.join(map(lambda x: str(float(x)), list(listi[:3])))
        for i in range(len('xyz')):                
            if listi[3 +i] < 0:
                files[2*i].write(entry)
                files[2*i].write('\t' + str(-(i+1))+'\n')
            elif listi[3+i] > 0:
                files[2*i + 1].write(entry)
                files[2*i + 1].write('\t' + str(i+1)+'\n')
        for ien, en in enumerate(energies):
            entry = [str(1e-3*en)+ '\t' + str(listi[3]) + '\t0\t0\t' + str(resi[ien]) + '\n',
                     str(1e-3*en)+ '\t0\t' + str(listi[3]) + '\t0\t' + str(resi[ien]) + '\n',
                     str(1e-3*en)+ '\t0\t0\t' + str(listi[3]) + '\t' + str(resi[ien]) + '\n']
            for i in range(len('xyz')):                
                if listi[3 +i] < 0:
                    files[2*i].write(entry[i])
                elif listi[3+i] > 0:
                    files[2*i + 1].write(entry[i])
    for dir_, dir_cnt in enumerate(directns.flatten()):
        if dir_cnt > 0:
            files[dir_].close()


def process_d_fl(name, energies_count, directns, lstf, totf, resf, varf):

    ptcl = 0
    files = {}
    for dir_, dir_cnt in enumerate(directns.flatten()):
        if dir_cnt > 0:
            files[dir_] = open(name + '_' + str(dir_) + '_' + str(ptcl) + ".spc", 'w')
            print name + str(dir_) + '_' + str(ptcl) + ".spc"

    resf.resize((len(resf)/energies_count), energies_count, 5)
    print resf[0]
    return
    for listi, toti, resi, vari in zip(lstf, totf, resf, varf):

        entry = '\t'.join(map(lambda x: str(float(x)), list(listi[:3])))
        for i in range(len('xyz')):
            if listi[3 + i] < 0:
                files[2 * i].write(entry)
                files[2 * i].write('\t' + str(-(i + 1)) + '\n')
            elif listi[3 + i] > 0:
                files[2 * i + 1].write(entry)
                files[2 * i + 1].write('\t' + str(i + 1) + '\n')
        for ien, en in enumerate(energies):
            entry = [str(1e-3 * en) + '\t' + str(listi[3]) + '\t0\t0\t' + str(resi[ien]) + '\n',
                     str(1e-3 * en) + '\t0\t' + str(listi[3]) + '\t0\t' + str(resi[ien]) + '\n',
                     str(1e-3 * en) + '\t0\t0\t' + str(listi[3]) + '\t' + str(resi[ien]) + '\n']
            for i in range(len('xyz')):
                if listi[3 + i] < 0:
                    files[2 * i].write(entry[i])
                elif listi[3 + i] > 0:
                    files[2 * i + 1].write(entry[i])
    for dir_, dir_cnt in enumerate(directns.flatten()):
        if dir_cnt > 0:
            files[dir_].close()

            
def process_det(path, config, d):
    lstf, totf, resf, varf = read_det_info(path, config, d)

    if 'FLUX' in d['type']:
        tmpl = read_tmpl(os.path.join(path, config["inidir"], d['templ']))
        if 'BASIC' in d['measure']:
            directns = count_directns_b(lstf)
            process_b_fl(d['name'], tmpl, directns, lstf, totf, resf, varf)
        elif 'DETAILED' in d['measure']:
            directns = count_directns_d(resf)
            process_b_fl(d['name'], len(tmpl['energies']), directns, lstf, totf, resf, varf)
    elif 'ENERGY' in d['type']:

        trees, energies = spc_en.filter_dets(grd, space, lays, totf, lstf)

        spc_en.process_en(grd, space, trees, en_space_distr, energies)
        
    elif 'TOK' in d['type']:
        trees, tok = spc_tok.filter_dets(grd, space, lays, totf, lstf)

        spc_tok.process_tok_x(grd, space, trees, jx_space_distr, tok[0])
        spc_tok.process_tok_y(grd, space, trees, jy_space_distr, tok[1])
        spc_tok.process_tok_z(grd, space, trees, jz_space_distr, tok[2])


def write_en_distr(path):
    with open(os.path.join(path, r'energy_distribution'), 'w') as f:
        for ix in range(len(grd['x']['i'])-1):
            for iy in range(len(grd['y']['i'])-1):
                for iz in range(len(grd['z']['i'])-1):
                    if abs(jx_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix+1, iy+1, iz+1]))+"\t{:.3E}\n".format(en_space_distr[ix, iy, iz]))


def write_jx_distr(path):
    with open(os.path.join(path, r'jx'), 'w') as f:
        for ix in range(len(grd['x']['i']) - 1):
            for iy in range(len(grd['y']['i'])):
                for iz in range(len(grd['z']['i'])):
                    if abs(jx_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix + 1, iy, iz])) + "\t{:.3E}\n".format(jx_space_distr[ix, iy, iz]))


def write_jy_distr(path):
    with open(os.path.join(path, r'jy'), 'w') as f:
        for ix in range(len(grd['x']['i'])):
            for iy in range(len(grd['y']['i']) - 1):
                for iz in range(len(grd['z']['i'])):
                    if abs(jy_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix, iy + 1, iz])) + "\t{:.3E}\n".format(jy_space_distr[ix, iy, iz]))


def write_jz_distr(path):
    with open(os.path.join(path, r'jz'), 'w') as f:
        for ix in range(len(grd['x']['i'])):
            for iy in range(len(grd['y']['i'])):
                for iz in range(len(grd['z']['i']) - 1):
                    #                    counter += 1
                    #                    if counter%20000 == 0:
                    #                        print '#',
                    if abs(jz_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix, iy, iz + 1])) + "\t{:.3E}\n".format(jz_space_distr[ix, iy, iz]))


if __name__ == '__main__':

    window = tk.Tk()
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = raw_input("Where's REMP?")
    # r'D:\ipm_prj\calculations\tzp_2\pechs'
    path = r'D:\ipm_prj\calculations\tzp_2'

    config, _, _, dets = read_PECHS.readPECHS(os.path.join(path, r'pechs'))
                       
    prj, grd_, space = read_REMP.readREMP(path)
    # triangles = geom(path)
    grd = {}
    grd['x'] = read_REMP.processGrid(grd_['x'])
    grd['y'] = read_REMP.processGrid(grd_['y'])
    grd['z'] = read_REMP.processGrid(grd_['z'])
    
    lays = np.unique(space.flatten())

    en_space_distr = np.full((len(grd['x']['i'])-1,
                              len(grd['y']['i'])-1,
                              len(grd['z']['i'])-1), 1e-32)

    jx_space_distr = np.full((len(grd['x']['i'])-1,
                              len(grd['y']['i']),
                              len(grd['z']['i'])), 0)
    jy_space_distr = np.full((len(grd['x']['i']),
                              len(grd['y']['i'])-1,
                              len(grd['z']['i'])), 0)
    jz_space_distr = np.full((len(grd['x']['i']),
                              len(grd['y']['i']),
                              len(grd['z']['i'])-1), 0)
    for det in dets:
        process_det(os.path.join(path, r'kuvsh2\kuvsh2'), config, det)

    write_en_distr(os.path.abspath(os.path.join(path, '..')), )
    write_jx_distr(os.path.abspath(os.path.join(path, '..')), )
    write_jy_distr(os.path.abspath(os.path.join(path, '..')), )
    write_jz_distr(os.path.abspath(os.path.join(path, '..')), )

