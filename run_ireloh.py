import importlib as imp
import os,glob
from crystals import Crystal
from utils import*
from import_cif import show_grid
import dutil                                ; imp.reload(dutil)
import multislice.mupy_utils as mut         ; imp.reload(mut)
import multislice.postprocess as pp         ; imp.reload(pp)
import multislice.multislice as mupy        ; imp.reload(mupy)
plt.close('all')

opath = '/home/tarik/Documents/git/ccp4/docs/figures/multislice/ireloh/ireloh2/'
# opath = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/ireloh/ireloh2_pad2/'
rpath = '/home/tarik/Documents/git/ccp4/reports/ireloh/xyz_2/'
path = '/home/tarik/Documents/git/ccp4/src/david_scripts/ED_processing/IRELOH_ED_Dataset_2-dials/'
datpath = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/ireloh/ireloh2_pad2/'
name = ''
step=1000

# dutil.gen_xyz(path,name=name,rep=20,pad=2,step=step,opath=datpath)
files = [os.path.basename(f) for f in glob.glob(datpath+'*.xyz')]
def run_simus():
    for i,data in enumerate(files):
        i_str = data.split('_')[-1].replace('.xyz','')  #;print(i_str)
        multi = mupy.Multislice(datpath,data=data,tail='%s' %i_str,cif_file='ireloh.cif',
            mulslice=False,keV=200,NxNy=2**12,repeat=[1,1,1],slice_thick=1,
            opt='srw',fopt='f',ssh='badb')

# run_simus()

def pp_simus():
    for i,data in enumerate(files):
        i_str = data.split('_')[-1].replace('.xyz','')
        pkl = datpath+'ireloh2_pad2_'+data.replace('.xyz','')+'_autoslic.pkl'
        multi = pp.load_multi_obj(pkl)
        multi.outf['patternS'] = multi.outf['patternnpy'].replace('.npy','S.npy')
        multi.postprocess(ppopt='uwS',ssh_alias='badb')
        qx,qy,I = np.load(multi._outf('patternS')); print(I.shape)
        # I/=I.max()
        # dsp.stddisp(im=[qx,qy,I],caxis=[0,1e-6],title=i_str,#,rings=[0.5,1,2],lw=2,
        #     xylims=0.3,imOpt='cv',axPos='V',cmap='binary',
        #     name=opath+'ireloh2_%ssim.png' %i_str,opt='sc')
        # # print(multi.check_simu_state())
pp_simus()

# dutil.show_phi(path)
# mut.show_grid(path+'ireloh2_0000.xyz','xy',opt='p')#xylims=[100,450,100,450])
# xylims={'xy':[0,450,0,450],'yz':[0,450,0,160],'xz':[0,450,0,160]}
# xylims=dict(zip(['xy','xz','yz'],[[0,450,0,450]]*3))
# xylims=dict(zip(['xy','xz','yz'],[500]*3))
# dutil.xyz_gif(path,rpath,opath,name=name,xylims=xylims)    #; p.wait()





def from_F():
    multi = mupy.Multislice(path,data='ireloh001.xyz',cif_file='ireloh.cif',
        mulslice=False,keV=200,NxNy=2**12,repeat=[1,1,1],slice_thick=1,
        opt='')
    hklM = 15
    (qxF,qyF,qzF),Fhkl = multi.get_structure_factor(hkl=None, hklMax=hklM, sym=1, v='')
    I_F = np.abs(Fhkl)**2 # perfect Bragg condition and flat ewald sphere
    sI = np.s_[:,:,hklM] # slice in the [0,0,1] direction
    I  = I_F[sI]
    I /= I.max()
    im=add_noise(I)
    im = broaden_peaks(im)
    np.save('dat/im001.npy',im)
    dsp.stddisp(im=[im],xylims=[0,516,0,516],cmap='binary',caxis=[0,0.02],imOpt='cv',axPos='V')

def add_noise(I):
    x,y = np.meshgrid(np.arange(-258,258),np.arange(-258,258))
    r = np.sqrt(x**2+y**2);r[r==0]=1
    im = np.random.rand(516,516)/(5*r)
    im[18:-16:16,18:-16:16] = I
    return im

def broaden_peaks(im):
    x0,y0 = np.meshgrid(range(-7,8),range(-7,8))
    Pb = np.exp(-0.5*(x0**2+y0**2))
    i = range(18,516-16,16)
    i,j = np.meshgrid(i,i)
    for i0,j0 in zip(i.flatten(),j.flatten()):
        im0 = im[i0,j0]
        im[i0+x0,j0+y0] += im[i0,j0]*Pb
        im[i0,j0] = im0
        # im[i0+x,j0+y] = np.maximum(im0[i0+x,j0+y],tol*1e-2)
    return im
# from_F()
