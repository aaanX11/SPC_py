import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.spatial import KDTree

#cells = np.loadtxt('E:\\ks_work\\calc\\kuvsh_emk_1\\KUVSH.CEL', np.dtype(int))#
cells = np.loadtxt('D:\\ipm_prj\\calculations\\kuvsh\\KUVSH.CEL', np.dtype(int))

nx = 116
ny = 116
nz = 200
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

gridF = open('D:\\ipm_prj\\calculations\\kuvsh\\KUVSH.GRD', 'r')
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





space.resize((nx,ny, nz))
print np.squeeze(space[:,ny/2:ny/2+1:1,:], axis=1).shape
print len(xi), len(yi),len(zi)

def get_range():
    xmin = 100
    xmax = -100
    ymin = 100
    ymax = -100
    zmin = 100
    zmax = -100
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if space[ix,iy,iz] ==4:
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
    enDistr = np.loadtxt('D:\\ipm_prj\\calculations\\emk_5\\energy_distribution')
    enSpaceCheck = np.full((nx, ny, nz), 1e-32)
    for entry in enDistr:
        ix, iy, iz = map(lambda x: int(x)-1,entry[:-1])
        enSpaceCheck[ix, iy, iz] = entry[-1]
    return enSpaceCheck
enSpaceCheck = f()
print np.unique(enSpaceCheck[:,58,:].T, axis = 0)
        
x05 = [0.5*(x1+x2) for x1,x2 in zip(xi[1:], xi[:-1])]
y05 = [0.5*(x1+x2) for x1,x2 in zip(yi[1:], yi[:-1])]
z05 = [0.5*(x1+x2) for x1,x2 in zip(zi[1:], zi[:-1])]
slicey = 58
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
im = ax.contourf(z05, x05,np.squeeze(enSpaceCheck[:,slicey,:]))
plt.colorbar(im)
plt.show()

enSpace = f()
slicez = 77
slicey = 76
slicex = 52
np.savetxt('xy.txt', np.squeeze(enSpace[xmin:xmax+1,ymin:ymax+1,slicez:slicez+1:1]), fmt='%1.2e')
np.savetxt('xz.txt', np.squeeze(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]), fmt='%1.2e')
print np.argmax(enSpace)
if np.isnan(enSpace).any():
    print 'alarm'
if np.isinf(enSpace).any():
    print 'aalarm'

#ind = np.unravel_index(np.argmax(a, axis=None), a.shape)
#https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html#numpy.argmax
fig = plt.figure()
ax = fig.add_subplot(2, 2, 1)
#slicez = nz/2
slicez = 77
im = ax.contourf(x05[xmin:xmax+1], y05[ymin:ymax+1],np.squeeze(enSpace[xmin:xmax+1,ymin:ymax+1,slicez:slicez+1:1]),
                 norm=colors.LogNorm(vmin=np.amin(enSpace[:,:,nz/2:nz/2+1:1]), vmax=np.amax(enSpace[:,:,slicez:slicez+1:1])),
                 cmap=plt.cm.coolwarm)
print np.amin(enSpace[:,:,nz/2:nz/2+1:1]), np.amax(enSpace[:,:,slicez:slicez+1:1])
plt.colorbar(im)
ax.set_title('energy distrib XY')
def addSubplot():
    ax = fig.add_subplot(2, 2, 2)
    #slicez = nz/2
    slicez = 77
    im = ax.contourf(x05[xmin:xmax+1], y05[ymin:ymax+1],np.squeeze(enSpaceCheck[xmin:xmax+1,ymin:ymax+1,slicez:slicez+1:1]),
                 norm=colors.LogNorm(vmin=np.amin(enSpaceCheck[:,:,nz/2:nz/2+1:1]), vmax=np.amax(enSpaceCheck[:,:,slicez:slicez+1:1])),
                 cmap=plt.cm.coolwarm)
    plt.colorbar(im)
    ax.set_title('energy distrib Check XY')
addSubplot()

def addAnotherSub():
    ax = fig.add_subplot(2, 2, 2)
    slicey = 76
    im = ax.contourf(z05[zmin:zmax+1], x05[xmin:xmax+1],np.squeeze(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]),
                     norm=colors.LogNorm(vmin=np.amin(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]), vmax=np.amax(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1])),                    
                     cmap=plt.cm.coolwarm)
    print np.amin(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1]), np.amax(enSpace[xmin:xmax+1,slicey:slicey+1:1,zmin:zmax+1])
    plt.colorbar(im)
    ax.set_title('energy distrib XZ')

addAnotherSub()
plt.show()
ax = fig.add_subplot(2, 2, 3)
slicex = 52
im = ax.contourf(z05[zmin:zmax+1], y05[ymin:ymax+1],np.squeeze(enSpace[slicex:slicex+1:1,ymin:ymax+1,zmin:zmax+1]),
                 norm=colors.LogNorm(vmin=np.amin(enSpace[slicex:slicex+1:1,ymin:ymax+1,zmin:zmax+1]), vmax=np.amax(enSpace[slicex:slicex+1:1,ymin:ymax+1,zmin:zmax+1])),                    
                 cmap=plt.cm.coolwarm)
print np.amin(enSpace[slicex:slicex+1:1,ymin:ymax+1,zmin:zmax+1]), np.amax(enSpace[slicex:slicex+1:1,ymin:ymax+1,zmin:zmax+1])
plt.colorbar(im)
ax.set_title('energy distrib YZ')

from mpl_toolkits.mplot3d import Axes3D
ax = fig.add_subplot(2, 2, 4, projection='3d')
im = ax.scatter(xx, yy, zz, c = enn)
plt.colorbar(im)

ax.set_title('energy distrib 3d')
plt.show()
#-1.5873535e+01   3.1478027e+01   2.0901855e+01 
