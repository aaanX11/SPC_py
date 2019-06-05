import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.spatial import KDTree
import os

def read_det_info(path, list_fn, tot_fn):
    lstf = np.loadtxt(list_fn)
    totf = np.loadtxt(tot_fn)

def read_energy_distribution(path):
    enDistr = np.loadtxt(os.path.join(path, r'energy_distribution'))
    enSpaceCheck = np.full((nx, ny, nz), 1e-32)
    for entry in enDistr:
        ix, iy, iz = map(lambda x: int(x)-1,entry[:-1])
        enSpaceCheck[ix, iy, iz] = entry[-1]
    return enSpaceCheck        


def filter_dets(grd, space, lays, totDet, listDet):
    
    counter = 0
    listDet_filtered = {l: [] for l in lays}
    totDet_filtered = {l: [] for l in lays}
    for tok, entry in zip(totDet, listDet[:, :-1]):
        if not (grd['x']['i'][0] <= entry[0] <= grd['x']['i'][-1] and grd['y']['i'][0] <= entry[1] <= grd['y']['i'][-1] and grd['z']['i'][0] <= entry[2] <= grd['z']['i'][-1]):
            continue
        
        ix1 = next((idx for idx, v in enumerate(grd['x']['i']) if v > entry[0]), len(grd['x']['i']))
        ix = [ix1]
        if entry[0] in grd['x']['i']:
            ix.append(ix1-1)
            
        iy1 = next((idx for idx, v in enumerate(grd['y']['i']) if v > entry[1]), len(grd['y']['i']))
        iy = [iy1]
        if entry[1] in grd['y']['i']:
            iy.append(iy1-1)
            
        iz1 = next((idx for idx, v in enumerate(grd['z']['i']) if v > entry[2]), len(grd['z']['i']))
        iz = [iz1]
        if entry[2] in grd['z']['i']:
            iz.append(iz1-1)        
        lays_near = np.unique(space(np.ix_(ix, iy, iz))).flatten()
        for lay in lays_near:
            listDet_filtered[lay].append(entry)
            totDet_filtered[lay].append(tok)
    trees = {lay: KDTree(listDet_filtered[lay]) for lay in lays}
    return trees, totDet_filtered


def process_tok_x(grd, space, trees, jx_space_distr, totDet):
    
    counter = 0
    for ix, x in enumerate(grd['x']['i05'][1:-1]):
        for iy, y in enumerate(grd['y']['i']):
            for iz, z in enumerate(grd['z']['i']):

                counter += 1
                if counter % 20000 == 0:
                    print '#',

                lays_near = (space[ix+1][iy][iz], space[ix+1][iy+1][iz], space[ix+1][iy][iz+1], space[ix+1][iy+1][iz+1])
                lays_near = np.delete(lays_near, np.where(lays_near == -1))
                lays_counted = np.unique(lays_near, return_counts=True)
                lay = lays_counted[0][[np.argmax(lays_counted[1])]]

                d, i = trees[lay].query((x, y, z))

                jx_space_distr[ix+1, iy, iz] += totDet[lay][i]


def process_tok_y(grd, space, trees, jy_space_distr, totDet):
    counter = 0
    for ix, x in enumerate(grd['x']['i']):
        for iy, y in enumerate(grd['y']['i05'][1:-1]):
            for iz, z in enumerate(grd['z']['i']):

                counter += 1
                if counter % 20000 == 0:
                    print '#',

                lays_near = (space[ix][iy + 1][iz], space[ix + 1][iy + 1][iz],
                             space[ix][iy + 1][iz + 1], space[ix + 1][iy + 1][iz + 1])
                lays_near = np.delete(lays_near, np.where(lays_near == -1))
                lays_counted = np.unique(lays_near, return_counts=True)
                lay = lays_counted[0][[np.argmax(lays_counted[1])]]

                d, i = trees[lay].query((x, y, z))

                jy_space_distr[ix, iy + 1, iz] += totDet[lay][i]


def process_tok_z(grd, space, trees, jz_space_distr, totDet):
    counter = 0
    for ix, x in enumerate(grd['x']['i']):
        for iy, y in enumerate(grd['y']['i']):
            for iz, z in enumerate(grd['z']['i05'][1:-1]):

                counter += 1
                if counter % 20000 == 0:
                    print '#',

                lays_near = (space[ix][iy][iz + 1], space[ix + 1][iy][iz + 1],
                             space[ix][iy + 1][iz + 1], space[ix + 1][iy + 1][iz + 1])
                lays_near = np.delete(lays_near, np.where(lays_near == -1))
                lays_counted = np.unique(lays_near, return_counts=True)
                lay = lays_counted[0][[np.argmax(lays_counted[1])]]

                d, i = trees[lay].query((x, y, z))

                jz_space_distr[ix, iy, iz + 1] += totDet[lay][i]


if __name__ == '__main__':

    enSpace = np.full((len(grd['x']['i'])-1, len(grd['y']['i'])-1, len(grd['z']['i']))-1, 1e-32)
    trees, totDet = filter_dets()
    get_this_to_work(enSpace)    

    slicez = len(grd['z']['i'])/2
    slicey = len(grd['y']['i'])/2
    slicex = len(grd['x']['i'])/2

    print np.argmax(enSpace)
    if np.isnan(enSpace).any():
        print 'alarm: NaN'
    if np.isinf(enSpace).any():
        print 'alarm: infty'

    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 1)

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

