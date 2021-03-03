import os
import numpy as np
from crystals import Crystal
from colorama import Fore as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
import utils.displayStandards as dsp

def import_cif(cif_file,xyz_file=None,rep=[1,1,1],pad=0):
    ''' convert cif file into autoslic .xyz input file
    - cif_file : file to import
    - rep : super cell repeat
    - pad : amount of padding on each side (in unit of super cell size)
    returns :
    - saves cif_file.xyz
    '''
    crys = Crystal.from_cif(cif_file)
    lat_params = np.array(crys.lattice_parameters[:3])
    Za     = np.array([a.atomic_number for a in crys.atoms])
    occ    = np.array([a.occupancy for a in crys.atoms])
    bfact  = np.ones(Za.shape)
    coords = np.array([a.coords_cartesian for a in crys.atoms])
    Nx,Ny,Nz = rep
    ax,by,cz = lat_params

    #replicate
    if sum(rep)>3 :
        ni,nj,nk = np.meshgrid(range(Nx),range(Ny),range(Nz))
        ni,nj,nk = ni.flatten(),nj.flatten(),nk.flatten()
        a1,a2,a3 = np.array(crys.lattice_vectors)
        coords = np.vstack([coords + i*a1+j*a2+k*a3 for i,j,k in zip(ni,nj,nk)])
        Za     = np.tile(Za   ,[ni.size])
        occ    = np.tile(occ  ,[ni.size])
        bfact  = np.tile(bfact,[ni.size])
        lat_params[0] = Nx*ax
        lat_params[1] = Ny*by
        lat_params[2] = Nz*cz

    #apply padding
    if pad:
        coords[:,0] += Nx*ax*pad
        coords[:,1] += Ny*by*pad
        coords[:,2] += Nz*cz*pad
        lat_params[0] *= 1+2*pad
        lat_params[1] *= 1+2*pad
        lat_params[2] *= 1+2*pad

    #save
    # print(Za.shape,coords.shape)
    pattern = np.hstack([Za[:,None],coords,occ[:,None],bfact[:,None]])
    if not xyz_file:xyz_file = cif_file.replace('.cif','.xyz')
    make_xyz(xyz_file,pattern,lat_params,fmt='%.4f')



def make_xyz(xyz_file,pattern,lat_params,fmt='%.4f'):
    '''Creates the.xyz file from a given compound and orientation
    - name    : Full path to the file to save
    - pattern : Nx6 ndarray - Z,x,y,z,occ,wobble format
    - lat_params :lattice parameters ax,by,cz
    '''
    compound = os.path.basename(xyz_file)
    ax,by,cz = lat_params

    header = 'one unit cell of %s\n' %(compound)
    header+= ' '.join([fmt]*3) %(ax,by,cz)
    np.savetxt(xyz_file,pattern,footer='-1',header=header,fmt='%d '+' '.join([fmt]*5),comments='')
    print(colors.GREEN+"coords file saved : \n"+colors.YELLOW+xyz_file+colors.RESET)

def import_xyz(xyz_file):
    '''import .xyz file
    - xyz_file : file to import
    returns :
    - pattern    : [Za,coords,occ,bfactor]
    - lat_params : [ax,by,cz] for the supercell
    '''
    with open(xyz_file,'r') as f:l=list(map(lambda s:s.strip().split(' '),f.readlines()))
    lat_params = np.array(l[1],dtype=float)
    pattern = np.array(l[2:-1],dtype=float)
    return pattern,lat_params

def show_grid(xyz_file,opts='',ms=1,**kwargs):
    '''display an .xyz file content in a 2D plane
    xyz_file : .xyz file
    opt : str format for 2 plane 'x1x2' - 'xy','xz','yz','zx',...
    '''
    pattern,lat_params = import_xyz(xyz_file)
    cs = {1:tuple([0.75]*3),6:tuple([0.25]*3),7:(0,0,1),8:(1,0,0),16:(1,1,0),17:(0,1,1)}
    Z = np.array(pattern[:,0],dtype=int)
    C = [cs[E] for E in Z]


    if opts =='3d':
        fig = plt.figure()
        ax = plt.subplot(111,projection='3d')
        ax.scatter3D(pattern[:,1],pattern[:,2],pattern[:,3],s=5,c=C)
        set_axes_equal(ax)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        fig.show()

    elif 'x' in opts or 'y' in opts or 'z' in opts:
        xij = {'x':1,'y':2,'z':3}
        x1,x2 = opts
        i,j = xij[x1],xij[x2]

        # fig,ax = plt.subplots()
        # plt.scatter(pattern[:,i],pattern[:,j],ms,C)
        # ax.add_patch(Rectangle((0,0),lat_params[i-1],lat_params[j-1],linewidth=2,edgecolor='b',alpha=0.1))
        # ax.set_xlabel('$%s$'%x1,fontsize=20)
        # ax.set_ylabel('$%s$'%x2,fontsize=20)
        # fig.show()
        scat = [pattern[:,i],pattern[:,j],ms,C]
        patches=[Rectangle((0,0),lat_params[i-1],lat_params[j-1],linewidth=2,edgecolor='b',alpha=0.1)]
        dsp.stddisp(scat=scat,patches=patches,labs=['$%s$'%x1,'$%s$'%x2],**kwargs)
    else:
        print('opts should be xy, xz, yz, or 3d')

def set_axes_equal(ax):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

##################################################################################################
#### test
##################################################################################################
if __name__ == '__main__':
    plt.close('all')
    figpath = '/home/tarik/Documents/git/ccp4/docs/projects/multislice/figures/ireloh/'
    # cif_file = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/ireloh/ireloh.cif'
    cif_file = 'ireloh.cif'
    xyz_file = 'ireloh111.xyz'
    # xyz_file = cif_file.replace('.cif','.xyz')
    import_cif(cif_file,xyz_file=xyz_file,rep=[1,1,1],pad=0)
    # import_cif(cif_file,rep=[15,15,15],pad=0)
    # show_grid(xyz_file,opts='xy')
    # show_grid('ireloh_rotated484.xyz',opts='3d')
    # show_grid('ireloh_rotated484.xyz',opts='xz',ms=1, opt='sc',name=figpath+'484_X.png')
