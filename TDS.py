'''post processing simulated data for ireloh'''
import importlib as imp;
from utils import*
from multislice import multislice as mupy; imp.reload(mupy)
from multislice import mupy_utils as mut ; imp.reload(mut)
from multislice import postprocess as pp ; imp.reload(pp)
import dutil                             ; imp.reload(dutil)
plt.close('all')

figpath = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/docs_fig/ireloh/'
datpath = 'simus/ireloh2/TDS/'
exp_path = 'IRELOH_ED_Dataset_2-dials/'
# df_name = datpath+'df.pkl'
image=149
bfacts = {1:0.01,6:0.1,7:0.1,8:0.1}

# dutil.gen_xyz(exp_path, name=datpath+'ireloh', image=image,rep=15, pad=2,popt=0)
# dutil.gen_xyz(exp_path, name=datpath+'ireloh_bf2_', image=image,rep=15, pad=2,bfacts=bfacts,popt=0)
def run_simu():
    multi = mupy.Multislice(datpath,data='ireloh_bf2_0149.xyz',
        #tail='%s' %i_str,cif_file='ireloh.cif',
        TDS=True,n_TDS=4,T=300,tail='bf2',
        mulslice=False,keV=200,NxNy=2**10,repeat=[1,1,1],slice_thick=1,
        opt='srw',fopt='f',ssh='badb')

# run_simu()

pkl = datpath+'TDS_bf2_autoslic.pkl'
multi = pp.load_multi_obj(pkl)
multi.postprocess(ppopt='uwS',ssh_alias='badb')


# dutil.plot_npy(multi._outf('patternS'),log=True,
#     caxis=[-9,-6],xylims=0.2,
#     name=figpath+'TDS_log.png',opt='p')

dutil.plot_npy(multi._outf('patternS'),log=False,
    caxis=[0,1e-7],xylims=0.2,title='TDS',
    name=figpath+'TDS_lin.png',opt='p')
