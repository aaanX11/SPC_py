# -*- coding: cp1251 -*-

import os
import sys
import numpy as np

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize

import Tkinter as tk
import tkMessageBox
import ttk


from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import matplotlib.collections as clt
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D

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
                 vec = map(float, v.split())
                 d["energies"] = [0.5*(e1+ e2) for e1, e2 in zip(vec[1:], vec[:-1])]
    if 'energies' not in d:
        interv = (d['maxen'] - d['minen']) / d['res']
        d['energies'] = [d['minen'] + 0.5 * interv + i * interv for i in range(d['res'])]
    return d
    

def read_det_info(path, config, d):
    lstf = np.loadtxt(os.path.join(path, config["inidir"], d['list']), ndmin=2)
    totf = np.loadtxt(os.path.join(path, config["resdir"], d['name']+'.tot'), ndmin=1)
    resf = np.genfromtxt(os.path.join(path, config["resdir"], d['name']+'.res'))
    print resf.ndim
    if resf.ndim == 1:
        resf_ = np.zeros((1, resf.shape[0]))
        resf_[0,:] = resf.copy()
        
        resf = resf_

    varf = np.loadtxt(os.path.join(path, config["resdir"], d['name']+'.var'), ndmin=1)
    return lstf, totf, resf, varf


def count_directns_b(lstf):
    directns = np.zeros((3, 2), dtype=int)
    
    for lst in lstf:
        if lst[3] > 0:
            directns[0, 1] += 1
        elif lst[3] < 0:
            directns[0, 0] += 1
        if lst[4] > 0:
            directns[1, 1] += 1
        elif lst[4] < 0:
            directns[1, 0] += 1
        if lst[5] > 0:
            directns[2, 1] += 1
        elif lst[5] < 0:
            directns[2, 0] += 1
    return directns


def count_directns_d(resf):
    directns = np.zeros((3, 2), dtype=int)
    maskx = (abs(resf[:, 1]) > 1e-32)
    # = np.unique(maskx, return_counts = True)
    _, ind = np.unique(maskx, return_counts=True)
    directns[0, :] = ind
    masky = (abs(resf[:, 2]) > 1e-32)
    _, ind = np.unique(maskx, return_counts=True)
    directns[1, :] = ind
    print directns[1, :]
    maskz = (abs(resf[:, 3]) > 1e-32)
    _, ind = np.unique(maskz, return_counts=True)
    directns[2, :] = ind
    print directns[2, :]
    return directns


def cap(f, particle, spc_index, en_len, dir_cnt):
    ptcls = {'ELECTRON': u'Электроны\n', 'POSITRON': u'Позитроны\n', 'PHOTON': u'Фотоны\n',
             'ELECTRONS': u'Электроны\n', 'POSITRONS': u'Позитроны\n', 'PHOTONS': u'Фотоны\n'}
    
    for p in ptcls:
        if p in particle.upper():
            particle_norm = p
    f.write(ptcls[particle_norm].encode('cp1251'))
    f.write(u'Номер спектра\n'.encode('cp1251'))
    f.write(str(spc_index))
    f.write('\n')
    f.write(u"Мощность спектра (шт/см**2/с)- для поверхностных, (шт/см**3/с)- для объемных\n1".encode('cp1251'))
    f.write(u"\nТип спектра (0-фиксированный, 1-разыгрывание, 2-список, 3-детекторы)\n3\nЧисло частиц\n".encode('cp1251'))
    f.write(str(en_len))
    f.write('\n')
    f.write(u"Количество детекторов\n".encode('cp1251'))
    f.write(str(dir_cnt))
    f.write('\n')
    f.write(u"Энергия+нормаль\n".encode('cp1251'))


def process_b_fl(name, tmpl, directns, lstf, totf, resf, varf):
    energies = tmpl['energies']

    ptcl = 0
    files = {}
    for dir_, dir_cnt in enumerate(directns.flatten()):
        #if dir_cnt > 0:
        files[dir_] = open(name + '_' + str(dir_) + '_' + str(ptcl) + ".spc", 'w')
        particle = name.split('_')[0]
        spc_ind = map(int, name.split('_')[1:])
        spc_ind += [dir_, ptcl]
        spc_ind = reduce(lambda a, x: 100 * a + x, spc_ind)
        cap(files[dir_], particle, spc_ind, len(energies), dir_cnt)
        print name + '_' + str(dir_) + '_' + str(ptcl) + ".spc"

    for listi, toti, resi, vari in zip(lstf, totf, resf, varf):
        if len(list(resi)) < len(energies):
            assert('Number of entries in detector is less than expected. Check detector type')
            
        entry = '\t'.join(map(lambda x: str(float(x)), list(listi[:3])))
        for i in range(len('xyz')):                
            if listi[3 + i] < 0:
                files[2*i].write(entry)
                files[2*i].write('\t' + str(-(i+1))+'\n')
            elif listi[3 + i] > 0:
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
            files[dir_] = open(os.path.join(path, name + '_' + str(dir_) + '_' + str(ptcl) + ".spc"), 'w')
            particle = name.split('_')[0]
            spc_ind = map(int, name.split('_')[1:])
            spc_ind += [dir_, ptcl]
            spc_ind = reduce(lambda a, x: 100 * a + x, spc_ind)
            print 'spc_ind = ', spc_ind
            cap(files[dir_], particle, spc_ind, energies_count, dir_cnt)
            print name + '_' + str(dir_) + '_' + str(ptcl) + ".spc"

    resf.resize((len(resf)/energies_count), energies_count, 5)

    for listi, toti, resi, vari in zip(lstf, totf, resf, varf):

        entry = '\t'.join(map(lambda x: str(float(x)), list(listi[:3])))
        for i in range(len('xyz')):
            if listi[3 + i] < 0:
                files[2 * i].write(entry)
                files[2 * i].write('\t' + str(-(i + 1)) + '\n')
            elif listi[3 + i] > 0:
                files[2 * i + 1].write(entry)
                files[2 * i + 1].write('\t' + str(i + 1) + '\n')
        for item in resi:
            entry = [str(1e-3 * item[0]) + '\t' + str(item[1]) + '\t0\t0\t' + str(item[4]) + '\n',
                     str(1e-3 * item[0]) + '\t0\t' + str(item[2]) + '\t0\t' + str(item[4]) + '\n',
                     str(1e-3 * item[0]) + '\t0\t0\t' + str(item[3]) + '\t' + str(item[4]) + '\n']
            for i in range(len('xyz')):
                if listi[3 + i] < 0:
                    files[2 * i].write(entry[i])
                elif listi[3 + i] > 0:
                    files[2 * i + 1].write(entry[i])
    for dir_, dir_cnt in enumerate(directns.flatten()):
        if dir_cnt > 0:
            files[dir_].close()


def get_energy_scatter(totf, lstf, varf):
    xx_ = []
    yy_ = []
    zz_ = []
    enn = []
    for en, entry, var in zip(totf, lstf, varf):
        if en < 1e-15 or var > 0.5:
            continue

        x_, y_, z_ = map(float, entry[:-1])
        xx_.append(x_)
        yy_.append(y_)
        zz_.append(z_)
        enn.append(en)

    return xx_, yy_, zz_, enn


def get_current_scatter(totf, lstf, varf):
    xx_ = []
    yy_ = []
    zz_ = []
    tok_x = []
    tok_y = []
    tok_z = []
    ampl = []
    for tok, location, var in zip(totf, lstf, varf):
        if var > 0.5:
            continue
        x_, y_, z_ = map(float, location[:-1])
        jx, jy, jz = map(float, tok)

        j = np.sqrt(jx*jx + jy*jy + jz*jz)

        xx_.append(x_)
        yy_.append(y_)
        zz_.append(z_)
        tok_x.append(jx/j)
        tok_y.append(jx/j)
        tok_z.append(jx/j)
        ampl.append(j)

    return xx_, yy_, zz_, tok_x, tok_y, tok_z, ampl


def process_det(path, config, d):
    #global ene
    lstf, totf, resf, varf = read_det_info(path, config, d)

    print d
    if 'FLUX' in d['type']:
        tmpl = read_tmpl(os.path.join(path, config["inidir"], d['templ']))
        if 'BASIC' in d['measure']:
            print 'areyoulistening?'
            directns = count_directns_b(lstf)
            process_b_fl(d['name'], tmpl, directns, lstf, totf, resf, varf)
        elif 'DETAILED' in d['measure']:
            directns = count_directns_b(lstf)
            process_d_fl(d['name'], len(tmpl['energies']), directns, lstf, totf, resf, varf)
    elif 'ENERGY' in d['type']:
        global xx_energy, yy_energy, zz_energy, enn_energy
        trees, energies = spc_en.filter_dets(grd, space, lays, totf, lstf)

        spc_en.process_en(grd, space, trees, en_space_distr, energies)

        xx_, yy_, zz_, enn = get_energy_scatter(totf, lstf, varf)
        xx_energy.append(xx_)
        yy_energy.append(yy_)
        zz_energy.append(zz_)
        enn_energy.append(enn)
        
    elif 'CURRENT' in d['type']:
        global xx_tok, yy_tok, zz_tok, tok_x, tok_y, tok_z, tok_values
        print 'CURRENT'
        trees, tok = spc_tok.filter_dets(grd, space, lays, totf, lstf)

        spc_tok.process_tok_x(grd, space, trees, jx_space_distr, tok)
        spc_tok.process_tok_y(grd, space, trees, jy_space_distr, tok)
        spc_tok.process_tok_z(grd, space, trees, jz_space_distr, tok)

        xx_, yy_, zz_, tok_x_, tok_y_, tok_z_, values = get_current_scatter(totf, lstf, varf)

        xx_tok.append(xx_)
        yy_tok.append(yy_)
        zz_tok.append(zz_)
        tok_x.append(tok_x_)
        tok_y.append(tok_y_)
        tok_z.append(tok_z_)
        tok_values.append(values)


def write_en_distr(path):
    if not np.any((abs(en_space_distr) > 1.e-32).flatten()):
        return
    with open(os.path.join(path, r'energy_distribution'), 'w') as f:
        for ix in range(len(grd['x']['i'])-1):
            for iy in range(len(grd['y']['i'])-1):
                for iz in range(len(grd['z']['i'])-1):
                    if abs(en_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix+1, iy+1, iz+1]))+"\t{:.3E}\n".format(en_space_distr[ix, iy, iz]))


def write_jx_distr(path):
    if not np.any((abs(jx_space_distr) > 1.e-32).flatten()):
        return
    with open(os.path.join(path, r'JX'), 'w') as f:
        for ix in range(len(grd['x']['i']) - 1):
            for iy in range(len(grd['y']['i'])):
                for iz in range(len(grd['z']['i'])):
                    if abs(jx_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix + 1, iy, iz])) + "\t{:.3E}\n".format(jx_space_distr[ix, iy, iz]))


def write_jy_distr(path):
    if not np.any((abs(jy_space_distr) > 1.e-32).flatten()):
        return
    with open(os.path.join(path, r'JY'), 'w') as f:
        for ix in range(len(grd['x']['i'])):
            for iy in range(len(grd['y']['i']) - 1):
                for iz in range(len(grd['z']['i'])):
                    if abs(jy_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix, iy + 1, iz])) + "\t{:.3E}\n".format(jy_space_distr[ix, iy, iz]))


def write_jz_distr(path):
    if not np.any((abs(jz_space_distr) > 1.e-32).flatten()):
        print 'why?'
        return
    print os.path.join(path, r'JZ')
    with open(os.path.join(path, r'JZ'), 'w') as f:
        for ix in range(len(grd['x']['i'])):
            for iy in range(len(grd['y']['i'])):
                for iz in range(len(grd['z']['i']) - 1):
                    #                    counter += 1
                    #                    if counter%20000 == 0:
                    #                        print '#',
                    if abs(jz_space_distr[ix, iy, iz]) > 1e-32:
                        f.write("\t".join(map(str, [ix, iy, iz + 1])) + "\t{:.3E}\n".format(jz_space_distr[ix, iy, iz]))


class Gui():
    def __init__(self):
        self.key_axis = 'x'
        self.fig = plt.figure()

        vmin_ = np.amin(en_space_distr)
        vmax_ = np.amax(en_space_distr)

        ax1 = self.fig.add_subplot(221) # plt.subplot2grid((12, 12), (0, 6), rowspan=5, colspan=5, fig=self.fig)
        ax2 = self.fig.add_subplot(222)
        ax3 = self.fig.add_subplot(223)
        ax1.set_xlabel('Y')
        ax1.set_ylabel('Z')
        ax2.set_xlabel('Z')
        ax2.set_ylabel('X')
        ax3.set_xlabel('Y')
        ax3.set_ylabel('X')

        self.data = en_space_distr[:, :, :]

        self.index = {'x': len(grd['x']) / 2, 'y': len(grd['y']) / 2, 'z': len(grd['z']) / 2}
        print self.index

        xx, yy = np.meshgrid(grd['z']['i'], grd['y']['i'])
        debug = False
        if debug:
            ax4 = self.fig.add_subplot(224)
            ax4.pcolormesh(xx, yy, space[58, 1:-1, 1:-1])
        xx, yy = np.meshgrid(grd['z']['i05'][1:-1], grd['y']['i05'][1:-1])
        self.im1 = ax1.pcolormesh(xx, yy, self.data[self.index['x'], :, :],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')
        xx, yy = np.meshgrid(grd['z']['i05'][1:-1], grd['x']['i05'][1:-1])
        self.im2 = ax2.pcolormesh(xx, yy, self.data[:, self.index['y'], :],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')
        xx, yy = np.meshgrid(grd['y']['i05'][1:-1], grd['x']['i05'][1:-1])
        self.im3 = ax3.pcolormesh(xx, yy, self.data[:, :, self.index['z']],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')

        self.cmap = plt.cm.get_cmap('tab20', 9 + 1)

        plt.colorbar(self.im1)

        if not debug and len(xx_energy) > 0:
            ax4 = self.fig.add_subplot(2, 2, 4, projection='3d')
            det_idx = 0
            ss = ax4.scatter(xx_energy[det_idx], yy_energy[det_idx], zz_energy[det_idx],
                             c=enn_energy[det_idx])
            plt.colorbar(ss)

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        plt.tight_layout()

    def on_key_press(self, event):
        if event.key in 'xyz':
            self.key_axis = event.key

    def on_scroll(self, event):
        print ("%s %s" % (event.button, event.step))
        dd = -1
        if event.button == 'up':
            dd = 1
        self.index[self.key_axis] = np.clip(self.index[self.key_axis] + int(event.step),
                                            0, len(grd[self.key_axis]['i'][1:-1]))
        print self.index
        self.update()

    def update(self):
        #if self.axis == 1:
        #    tt_ =self.X[self.ind, :, :], self.kr_)
        #elif self.axis == 2:
        #    tt_ = self.X[:, self.ind, :], self.kr_)
        #elif self.axis == 3:
        #    tt_ = self.X[:, :, self.ind], self.kr_)
        ##        print ("full Energy = {0}".format(np.max(tt_)))

        self.im1.set_array(self.data[self.index['x'], :, :].ravel())
        self.im2.set_array(self.data[:, self.index['y'], :].ravel())
        self.im3.set_array(self.data[:, :, self.index['z']].ravel())
        self.fig.canvas.draw()


class GuiTok():
    def __init__(self):
        self.key_axis = 'x'
        self.fig = plt.figure()

        self.ax1 = self.fig.add_subplot(221)
        self.ax2 = self.fig.add_subplot(222)
        self.ax3 = self.fig.add_subplot(223)
        self.ax1.set_xlabel('Y')
        self.ax1.set_ylabel('Z')
        self.ax2.set_xlabel('Z')
        self.ax2.set_ylabel('X')
        self.ax3.set_xlabel('Y')
        self.ax3.set_ylabel('X')

        self.index = {'x': len(grd['x']) / 2, 'y': len(grd['y']) / 2, 'z': len(grd['z']) / 2}
        print self.index

        self.data = {}
        self.set_data('x')
        self.set_pics()
        xx, yy = np.meshgrid(grd['z']['i'], grd['y']['i'])
        debug = False
        if debug:
            ax4 = self.fig.add_subplot(224)
            ax4.pcolormesh(xx, yy, space[58, 1:-1, 1:-1])

        self.cmap = plt.cm.get_cmap('tab20', 9 + 1)

        plt.colorbar(self.im1)

        if not debug and len(xx_tok) > 0:
            ax4 = self.fig.add_subplot(2, 2, 4, projection='3d')

            det_idx = 0
            ss = ax4.quiver(xx_tok[det_idx], yy_tok[det_idx], zz_tok[det_idx],
                            tok_x[det_idx], tok_y[det_idx], tok_z[det_idx],
                            tok_values[det_idx])
            #plt.colorbar(ss)

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        plt.tight_layout()
    def set_data(self, key):
        if key == 'x':
            self.data['j'] = jx_space_distr[:, :, :]
            self.data['x'] = grd['x']['i05'][1:-1]
            self.data['y'] = grd['y']['i'][:]
            self.data['z'] = grd['z']['i'][:]
        elif key == 'y':
            self.data['j'] = jx_space_distr[:, :, :]
            self.data['x'] = grd['x']['i05'][1:-1]
            self.data['y'] = grd['y']['i'][:]
            self.data['z'] = grd['z']['i'][:]
        elif key == 'z':
            self.data['j'] = jx_space_distr[:, :, :]
            self.data['x'] = grd['x']['i05'][1:-1]
            self.data['y'] = grd['y']['i'][:]
            self.data['z'] = grd['z']['i'][:]

        self.index['x'] = len(self.data['x']) / 2
        self.index['y'] = len(self.data['y']) / 2
        self.index['z'] = len(self.data['z']) / 2


    def set_pics(self):
        vmin_ = np.amin(self.data['j'])
        vmax_ = np.amax(self.data['j'])
        xx, yy = np.meshgrid(self.data['z'], self.data['y'])
        self.im1 = self.ax1.pcolormesh(xx, yy, self.data['j'][self.index['x'], :, :],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')
        xx, yy = np.meshgrid(self.data['z'], self.data['x'])
        self.im2 = self.ax2.pcolormesh(xx, yy, self.data['j'][:, self.index['y'], :],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')
        xx, yy = np.meshgrid(self.data['y'], self.data['x'])
        self.im3 = self.ax3.pcolormesh(xx, yy, self.data['j'][:, :, self.index['z']],
                                  vmin=vmin_, vmax=vmax_, shading='gouraud')

    def on_key_press(self, event):
        if event.key in 'xyz':
            self.key_axis = event.key

    def on_scroll(self, event):
        print ("%s %s" % (event.button, event.step))
        dd = -1
        if event.button == 'up':
            dd = 1

        if event.inaxis in self.ax1:
            key_axis = 'x'
        if event.inaxis in self.ax2:
            key_axis = 'y'
        if event.inaxis in self.ax3:
            key_axis = 'z'

        self.index[key_axis] = np.clip(self.index[key_axis] + int(event.step),
                                       0, len(self.data[key_axis]))
        self.update()

    def update(self):
        #if self.axis == 1:
        #    tt_ =self.X[self.ind, :, :], self.kr_)
        #elif self.axis == 2:
        #    tt_ = self.X[:, self.ind, :], self.kr_)
        #elif self.axis == 3:
        #    tt_ = self.X[:, :, self.ind], self.kr_)
        ##        print ("full Energy = {0}".format(np.max(tt_)))

        self.im1.set_array(self.data['j'][self.index['x'], :, :].ravel())
        self.im2.set_array(self.data['j'][:, self.index['y'], :].ravel())
        self.im3.set_array(self.data['j'][:, :, self.index['z']].ravel())
        self.fig.canvas.draw()


if __name__ == '__main__':

    #window = tk.Tk()
    if len(sys.argv) > 2:
        path = sys.argv[1]
        path_p = os.path.abspath(os.path.join(path, sys.argv[2]))
    else:
        path = raw_input("Where's REMP?")
        path_p = os.path.abspath(os.path.join(path, raw_input("Where's PECHS?")))
    # r'D:\ipm_prj\calculations\tzp_2\pechs'
    #path = r'D:\ipm_prj\calculations\tzp_2'
    #path_p = os.path.abspath(os.path.join(path, r'pechs'))

    config, _, _, dets = read_PECHS.readPECHS(path_p)
                       
    prj, grd_, space = read_REMP.readREMP(path)

    # triangles = geom(path)
    grd = {}
    grd['x'] = read_REMP.process_grid(grd_['x'])
    grd['y'] = read_REMP.process_grid(grd_['y'])
    grd['z'] = read_REMP.process_grid(grd_['z'])
    
    lays = np.unique(space.flatten())

    en_space_distr = np.full((len(grd['x']['i'])-1,
                              len(grd['y']['i'])-1,
                              len(grd['z']['i'])-1), 1e-32)

    jx_space_distr = np.full((len(grd['x']['i'])-1,
                              len(grd['y']['i']),
                              len(grd['z']['i'])), 0, dtype=np.dtype(float))
    jy_space_distr = np.full((len(grd['x']['i']),
                              len(grd['y']['i'])-1,
                              len(grd['z']['i'])), 0, dtype=np.dtype(float))
    jz_space_distr = np.full((len(grd['x']['i']),
                              len(grd['y']['i']),
                              len(grd['z']['i'])-1), 0, dtype=np.dtype(float))
    print dets
    xx_energy = []
    yy_energy = []
    zz_energy = []
    enn_energy = []

    xx_tok = []
    yy_tok = []
    zz_tok = []
    tok_x = []
    tok_y = []
    tok_z = []
    tok_values = []

    for det in dets:

        process_det(path_p, config, det)
        lstf, totf, resf, varf = read_det_info(path_p, config, det)

    write_en_distr(os.path.abspath(path))

    write_jx_distr(os.path.abspath(path))
    write_jy_distr(os.path.abspath(path))
    write_jz_distr(os.path.abspath(path))

    g = Gui()
    plt.show()

    g = GuiTok()
    plt.show()
