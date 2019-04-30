import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.spatial import KDTree

cells = np.loadtxt(r'D:\ipm_prj\calculations\tzp_2\KUVSH.CEL', np.dtype(int))

listDet = np.loadtxt(r'D:\ipm_prj\calculations\tzp_2\pechs\initials\volume_5.lst')
totDet = np.loadtxt(r'D:\ipm_prj\calculations\tzp_2\pechs\results\energy_5.res')



gridF = open(r'D:\ipm_prj\calculations\emk_5\KUVSH.GRD', 'r')
gridF.readline()#<FULL>
gridF.readline()#1
gridF.readline()#<EQUAL (0-not equal 1-equal)>
gridF.readline()#0
gridF.readline()#X
gridF.readline()#116
xi = map(float,gridF.readline().strip().split())
gridF.readline()#Y
gridF.readline()#116
yi = map(float,gridF.readline().strip().split())
gridF.readline()#Z
gridF.readline()#
zi = map(float,gridF.readline().strip().split())
gridF.close()


nx = len(xi)-1
ny = len(yi)-1
nz = len(zi)-1
space = np.array([-1]*nx*ny*nz)

print np.sum(cells[1::2])
print cells.size
print cells


filled = 0
for i in range(cells.size/2):
    fill = cells[2*i+1]
    space[filled:filled+fill] = cells[2*i]
    filled += fill
print filled




space.resize((nx,ny, nz))
print np.squeeze(space[:,ny/2:ny/2+1:1,:], axis=1).shape
print len(xi), len(yi),len(zi)

def get_range():
    xmin = xi[-1]-1
    xmax = xi[0]
    ymin = yi[-1]-1
    ymax = yi[0]
    zmin = zi[-1]-1
    zmax = zi[0]
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if space[ix,iy,iz] == 5:
                    if ix >xmax:
                        xmax = ix
                    if ix < xmin:
                        xmin = ix
                    if iy >ymax:
                        ymax = iy
                    if iy < ymin:
                        ymin = iy
                    if iz >zmax:
                        zmax = iz
                    if iz < zmin:
                        zmin = iz
    return xmin, xmax, ymin, ymax, zmin, zmax
xmin, xmax, ymin, ymax, zmin, zmax = get_range()
print xmin, xmax, ymin, ymax, zmin, zmax


def f():
    enDistr = np.loadtxt(r'D:\ipm_prj\calculations\tzp_2\energy_distribution')
    enSpaceCheck = np.full((nx, ny, nz), 1e-32)
    for entry in enDistr:
        ix, iy, iz = map(lambda x: int(x)-1,entry[:-1])
        enSpaceCheck[ix, iy, iz] = entry[-1]
    return enSpaceCheck
#enSpaceCheck = f()
        

xx = []
yy = []
zz = []
enn = []
for en,entry in zip(totDet,listDet):
    x,y,z = map(float,entry[:-1])
    if en >1e-15:
        xx.append(x)
        yy.append(y)
        zz.append(z)
        enn.append(en)
x05 = [0.5*(x1+x2) for x1,x2 in zip(xi[1:], xi[:-1])]
y05 = [0.5*(x1+x2) for x1,x2 in zip(yi[1:], yi[:-1])]
z05 = [0.5*(x1+x2) for x1,x2 in zip(zi[1:], zi[:-1])]


def filter_dets():
    counter = 0
    listDet_filtered = []
    totDet_filtered = []
    for en,entry in zip(totDet, listDet[:, :-1]):
        #if x05[xmin]<=entry[0]<=x05[xmax] and y05[ymin]<=entry[1]<=y05[ymax] and z05[zmin]<=entry[2]<=z05[zmax]:
        listDet_filtered.append(entry)
        totDet_filtered.append(en)
    return listDet_filtered, totDet_filtered
listDet, totDet = filter_dets()

def get_this_to_work():
    unique, counts = np.unique(space, return_counts=True)
    
    enSpaceToFile = np.zeros((dict(zip(unique, counts))[4],4))
    tree = KDTree(listDet)
    counter = 0
    with open(r'D:\ipm_prj\calculations\tzp_2\energy_distribution', 'w') as f:
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if space[ix,iy,iz] == 5:
                        
                        d,i = tree.query((x05[ix],y05[iy],z05[iz]))
                        
                        f.write("\t".join(map(str,[ix+1, iy+1, iz+1]))+"\t{:.3E}\n".format(totDet[i]))
    
get_this_to_work()
def energy_distribution():    
    enSpace = np.full((nx, ny, nz), 1e-32)
    tree = KDTree(listDet)
    counter = 0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if space[ix,iy,iz] == 5:
                    counter += 1
                    if counter%20000 == 0:
                        print '#',
                    
                    d,i = tree.query((x05[ix],y05[iy],z05[iz]))
                    
                    enSpace[ix, iy, iz] = totDet[i]
    return enSpace
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if enSpace[ix,iy,iz] >0:
                    if ix >xmax1:
                        xmax1 = ix
                    if ix < xmin1:
                        xmin1 = ix
                    if iy >ymax1:
                        ymax1 = iy
                    if iy < ymin1:
                        ymin1 = iy
                    if iz >zmax1:
                        zmax1 = iz
                    if iz < zmin1:
                        zmin1 = iz
    return enSpace
enSpace = energy_distribution()
slicez = nz/2
slicey = ny/2
slicex = nz/2
np.savetxt('xy.txt', np.squeeze(enSpace[:,:,slicez:slicez+1:1]), fmt='%1.2e')
np.savetxt('xz.txt', np.squeeze(enSpace[:,slicey:slicey+1:1,:]), fmt='%1.2e')
print np.argmax(enSpace)
if np.isnan(enSpace).any():
    print 'alarm'
if np.isinf(enSpace).any():
    print 'aalarm'

#ind = np.unravel_index(np.argmax(a, axis=None), a.shape)
#https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html#numpy.argmax
fig = plt.figure()
ax = fig.add_subplot(2, 2, 1)
slicez = nz/2
print len(xi), len(yi), np.squeeze(enSpace[:,:,slicez:slicez+1:1]).shape
im = ax.pcolor(yi, xi, np.squeeze(enSpace[:,:,slicez:slicez+1:1]),
                 norm=colors.LogNorm(vmin=np.amin(enSpace[:,:,nz/2:nz/2+1:1]), vmax=np.amax(enSpace[:,:,slicez:slicez+1:1])),
                 cmap=plt.cm.coolwarm)
print np.amin(enSpace[:,:,nz/2:nz/2+1:1]), np.amax(enSpace[:,:,slicez:slicez+1:1])
plt.colorbar(im)
ax.set_title('energy distrib XY')
def addSubplot():
    ax = fig.add_subplot(2, 2, 2)
    slicez = nz/2
    im = ax.pcolor(xi, yi,np.squeeze(enSpace[:,:,nz/2:nz/2+1:1]),
                 norm=colors.LogNorm(vmin=np.amin(enSpace[:,:,nz/2:nz/2+1:1]), vmax=np.amax(enSpace[:,:,nz/2:nz/2+1:1])),
                 cmap=plt.cm.coolwarm)
    plt.colorbar(im)
    ax.set_title('energy distrib Check XY')
addSubplot()

def addAnotherSub():
    ax = fig.add_subplot(2, 2, 2)
    slicey = ny/2
    im = ax.pcolor(zi, xi,np.squeeze(enSpace[:,slicey:slicey+1:1,:]),
                     norm=colors.LogNorm(vmin=np.amin(enSpace[:,slicey:slicey+1:1,:]), vmax=np.amax(enSpace[:,slicey:slicey+1:1,:])),                    
                     cmap=plt.cm.coolwarm)
    plt.colorbar(im)
    ax.set_title('energy distrib XZ')

ax = fig.add_subplot(2, 2, 3)
slicex = nx/2
im = ax.pcolor(zi, yi,np.squeeze(enSpace[slicex:slicex+1:1,:,:]),
                 norm=colors.LogNorm(vmin=np.amin(enSpace[slicex:slicex+1:1,:,:]), vmax=np.amax(enSpace[slicex:slicex+1:1,:,:])),                    
                 cmap=plt.cm.coolwarm)
plt.colorbar(im)
ax.set_title('energy distrib YZ')

from mpl_toolkits.mplot3d import Axes3D
ax = fig.add_subplot(2, 2, 4, projection='3d')
im = ax.scatter(xx, yy, zz, c = enn)
plt.colorbar(im)

ax.set_title('energy distrib 3d')
plt.show()
#-1.5873535e+01   3.1478027e+01   2.0901855e+01 
