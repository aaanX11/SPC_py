import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import Tkinter as tk
import tkMessageBox
import ttk



import matplotlib
import matplotlib.pyplot as plt

from matplotlib.widgets import Button
from matplotlib.widgets import TextBox

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def readStart(startFN):
    d={}
    with open(startFN) as startf:
        d["prjf"] = startf.readline().strip()
        #d["npx"] = int(startf.readline().strip())
    return d
def readPrj(prjFN):
    d={}
    with open(prjFN) as prjf:
        
        line = prjf.readline().strip()
        while line:
            line = prjf.readline().strip()
            #print line.strip("><")
            if "Dxf name" in line:
                d["dxff"] = prjf.readline().strip()
            elif "Grd name" in line:
                d["grdf"]= prjf.readline().strip()
            elif "Cell name" in line:
                d["celf"]= prjf.readline().strip()
    print d
    return d

def spaceOf(celFN, shape):
    cells = np.loadtxt(celFN, np.dtype(int))
    size = reduce(lambda x,y: x*y, shape)
    space = np.array([-1]*size)

    assert np.sum(cells[1::2]) == size
   
    filled = 0
    for i in range(cells.size/2):
        fill = cells[2*i+1]
        space[filled:filled+fill] = cells[2*i]
        filled += fill
        
    space.resize(shape)
    return space
    
def readGrd(grdFN):
    d = {}
    with open(grdFN) as grdf:
        grdf.readline() #<FULL>
        grdf.readline() #1
        grdf.readline() #<EQUAL (0-not equal 1-equal)>
        grdf.readline() #0
        grdf.readline() #X
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
def processGrid(xi):
    
    d = {}
    
    xi05 = [0.5*(x1+x2) for x1,x2 in zip(xi[1:], xi[:-1])]

    dxi05 = [x2-x1 for x1,x2 in zip(xi[1:], xi[:-1])]
    dxi05.insert(0,dxi05[0])
    dxi05.append(dxi05[-1])
    
    xi05.insert(0, xi[0]-0.5*dxi05[0])
    xi05.append(xi05[-1]+0.5*dxi05[-1])
    d['i'] = xi
    d['i05'] = xi05
    d['di05'] = dxi05
    
    return d

def readREMP(path):
    if os.path.exists(os.path.join(path, "START_N")):
        start = readStart(os.path.join(path, "START_N"))
    else:
        try:
            start = readStart(os.path.join(path, "START"))
        except (IOError, ValueError) as e: 
            print repr(e)
    prj = readPrj(os.path.join(path, start["prjf"]))
    grd = readGrd(os.path.join(path, prj["grdf"]))
    space = spaceOf(os.path.join(path, prj["celf"]), (grd['nx'],grd['ny'],grd['nz']))
    return prj, grd, space

def changeSubplot(i, ax, cells,x, y, cmap):
    xx, yy = np.meshgrid(y,x)
    print cells.shape
    print len(x), len(y)
    ax.pcolor(xx,yy,cells, cmap=cmap)

    
def addSubplot(i,fig, cells, x,y, cmap, lays):
    ax = fig.add_subplot(2, 2, i)
    #lays = np.unique(cells.flatten())
    print cells.shape
    print len(x), len(y)
    xx, yy = np.meshgrid(y,x)

    
    #USE PCOLORMESH INSTEAD
    im = ax.pcolor(xx,yy,cells, cmap=cmap)



    #im.cmap.set_under('yellow')
    #im.cmap.set_over('cyan')
    
    return ax

def addAnotherSub(fig, cells, x,y,z, cmap, lays):
    ax = fig.add_subplot(2, 2, 2)
    #lays = np.unique(cells.flatten())
    im = ax.contourf(z, x,np.squeeze(cells[:,len(y)/2,:]), lays, extend = 'both', cmap=cmap)
    #im.cmap.set_under('yellow')
    #im.cmap.set_over('cyan')
    ax.set_title('cells XZ central')
    return im
def addYetAnotherSub(fig, cells, x,y,z, cmap, lays):
    ax = fig.add_subplot(2, 2, 3)
    
    im = ax.contourf(z, y,np.squeeze(cells[len(x)/2,:,:]), lays, extend = 'both', cmap=cmap)
    #im.cmap.set_under('yellow')
    #im.cmap.set_over('cyan')
    ax.set_title('cells XZ central')
    return im
    
def add3dSub(fig, cells, x,y,z, cmap, lays):
    ax = fig.add_subplot(2, 2, 4)
    slicey = 76
    im = ax.contourf(z05[zmin:zmax+1], x05[xmin:xmax+1],np.squeeze(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]),
                     norm=colors.LogNorm(vmin=np.amin(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]), vmax=np.amax(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1])),                    
                     cmap=plt.cm.coolwarm)
    print np.amin(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]), np.amax(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1])
    plt.colorbar()
    ax.legend()
    ax.set_title('energy distrib XZ')

def on_key_press(event):#, canvas, toolbar):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)


class gui:
    def __init__(self, master, space, lays, grd):
        self.master = master

        self.frame = tk.Frame(self.master)
        self.frame.pack(side = tk.TOP)
        self.fig = plt.figure(figsize=(5, 4))
        
        self.cmap = plt.cm.get_cmap('tab20', np.max(lays)+1)
        self.space = space
        self.grd = grd

        self.x = len(grd['x']['i'])/2
        self.y = len(grd['y']['i'])/2
        self.z = len(grd['z']['i'])/2
        
        self.ax1 = addSubplot(1, self.fig, np.squeeze(self.space[:,:,self.z]), self.grd['x']['i'], grd['y']['i'], self.cmap, lays)
        self.ax1.set_title('cells XY central')
        self.ax1.format_coord = lambda x,y : self.format_coord1(x,y)
        
        self.ax2 = addSubplot(2,self.fig, np.squeeze(self.space[:,self.y,:]), self.grd['x']['i'], grd['z']['i'], self.cmap, lays)
        self.ax2.set_title('cells XZ central')
        self.ax2.format_coord = lambda x,y : self.format_coord2(x,y)
        
        self.ax3 = addSubplot(3,self.fig, np.squeeze(self.space[self.x,:,:]), self.grd['y']['i'], grd['z']['i'], self.cmap, lays)
        self.ax3.set_title('cells YZ central')
        self.ax3.format_coord = lambda x,y : self.format_coord3(x,y)
        
        plt.figlegend([mpatches.Patch(color=self.cmap(b)) for b in lays.flatten()],
                   map(lambda x : str(int(x)),lays.flatten()), loc = (0.5, 0))

        self.fig.canvas.set_window_title("SECTION")

        #plt.subplots_adjust(bottom=0.2)

        #self.nb = ttk.Notebook(self.master)
        #self.nb.grid(row = 0, column = 0, columnspan=4, rowspan=4, pady=0,padx=0)

        #page1 = ttk.Frame(self.nb)
        #page2 = ttk.Frame(self.nb)
        
        tk.Label(self.frame, text=u"POINT").grid(row=1, column=1)

        
        
        #self.xStr = tk.StringVar(value = str(grd['x']['i05'][self.x]))
        #self.yStr = tk.StringVar(value = str(grd['y']['i05'][self.y]))
        #self.zStr = tk.StringVar(value = str(grd['z']['i05'][self.z]))

        self.xStr = tk.StringVar(value = str(self.x))
        self.yStr = tk.StringVar(value = str(self.y))
        self.zStr = tk.StringVar(value = str(self.z))
        
        xBox = tk.Entry(self.frame, textvariable=self.xStr)            
        xBox.bind('<Return>', self.changeCoefs)
        xBox.grid(row=2, column=1)
        yBox = tk.Entry(self.frame, textvariable=self.yStr)            
        yBox.bind('<Return>', self.changeCoefs)
        yBox.grid(row=2, column=2)
        zBox = tk.Entry(self.frame, textvariable=self.zStr)            
        zBox.bind('<Return>', self.changeCoefs)
        zBox.grid(row=2, column=3)

        #self.nb.add(page1, text = 'PECHS')
        #self.nb.add(page2, text = 'PECHS+')
        #self.nb.select(page2)
        
        fig_canvas_agg = FigureCanvasTkAgg(self.fig, master = self.master)
        fig_canvas_agg.draw()
        fig_canvas_agg.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


        toolbar = NavigationToolbar2TkAgg(fig_canvas_agg, self.master)
        toolbar.update()
        #fig_canvas_agg.get_tk_widget().grid(row=5, column =0, columnspan=4)
        fig_canvas_agg._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

        fig_canvas_agg.mpl_connect("key_press_event", on_key_press)
        #saveButton = tk.Button( self.master,text = 'Save', command = self.saveTF)
        #saveButton.grid(row=7, column=2)        

    def format_coord1(self, x, y):
        col = next(i[0] for i in enumerate(self.grd['y']['i']) if i[1] > x)
        row = next(i[0] for i in enumerate(self.grd['x']['i']) if i[1] > y)
        numcols = len(grd['y']['i'])
        numrows = len(grd['x']['i'])
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            return 'x=%d, y=%d, z=%d' % (row, col, self.z)
        
    def format_coord2(self, x, y):
        col = next(i[0] for i in enumerate(self.grd['z']['i']) if i[1] > x)
        row = next(i[0] for i in enumerate(self.grd['x']['i']) if i[1] > y)
        numcols = len(grd['z']['i'])
        numrows = len(grd['x']['i'])
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            return 'x=%d, y=%d, z=%d' % (row, self.y, col)
        else:
            print 'ALAAAARM'
        
    def format_coord3(self, x, y):
        col = next(i[0] for i in enumerate(self.grd['z']['i']) if i[1] > x)
        row = next(i[0] for i in enumerate(self.grd['y']['i']) if i[1] > y)
        numcols = len(grd['z']['i'])
        numrows = len(grd['y']['i'])
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            return 'x=%d, y=%d, z=%d' % (self.x, row, col)
        else:
            print 'ALAAAARM'
       
    def changeCoefs(self, event):
        if not isNum(self.xStr.get()):
            self.xStr.set(self.x)
        else:
            i = int(self.xStr.get())
            if 0 <= i < len(grd['x']['i']):
                self.x = i
            else:
                self.xStr.set(self.x)
            
        if not isNum(self.yStr.get()):
            self.yStr.set(self.y)
        else:
            i = int(self.yStr.get())
            if 0 <= i < len(grd['y']['i']):
                self.y = i
            else:
                self.yStr.set(self.y)

        if not isNum(self.zStr.get()):
            self.zStr.set(self.z)
        else:
            i = int(self.zStr.get())
            if 0 <= i < len(grd['z']['i']):
                self.z = i
            else:
                self.zStr.set(self.z)

        changeSubplot(1, self.ax1, np.squeeze(self.space[:,:,self.z]), self.grd['x']['i'], grd['y']['i'], self.cmap)
        self.ax1.set_title('cells XY z = {}'.format(self.z))
        
        changeSubplot(2,self.ax2, np.squeeze(self.space[:,self.y,:]), self.grd['x']['i'], grd['z']['i'], self.cmap)
        self.ax2.set_title('cells XZ y = {}'.format(self.y))
        
        changeSubplot(3,self.ax3, np.squeeze(self.space[self.x,:,:]), self.grd['y']['i'], grd['z']['i'], self.cmap)
        self.ax3.set_title('cells YZ x = {}'.format(self.x))
        
        self.fig.canvas.draw_idle()

                
if __name__ == '__main__':
    window = tk.Tk()
    if(len(sys.argv)>1):
        path = sys.argv[1]
    else:
        path = raw_input("Where's REMP?")
                
    prj, grd_, space = readREMP(path)
    grd = {}
    grd['x'] = processGrid(grd_['x'])
    grd['y'] = processGrid(grd_['y'])
    grd['z'] = processGrid(grd_['z'])
    
    lays = np.unique(space.flatten())
    

    window.wm_title("SECTION")
    g = gui(window, space, lays, grd)
    window.protocol("WM_DELETE_WINDOW", lambda : window.quit())
    window.mainloop()
    window.destroy()   
