'''post processing simulated data for ireloh'''
import importlib as imp;
from utils import*
from multislice import mupy_utils as mut ; imp.reload(mut)
from multislice import postprocess as pp ; imp.reload(pp)
import dutil                             ; imp.reload(dutil)
plt.close('all')

figpath = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/docs_fig/ireloh/'
datpath = 'simus/ireloh2/pad_tests/'
exp_path = 'IRELOH_ED_Dataset_2-dials/'
df_name = datpath+'df.pkl'
image=149

pads = [1,2,3]*3
reps = [10]*3 +[15]*3+[20]*3

def gen_simus():
    for rep,pad in zip(reps,pads):
        dutil.gen_xyz(exp_path, name=datpath+'pad%d-rep%d_'%(pad,rep), image=image,
            rep=rep, pad=pad, popt=0)

def run_simus():
    vals = ['pad%d-rep%d_%s.xyz'%(pad,rep,str(image).zfill(4)) for pad,rep in zip(pads,reps)]
    df = mut.sweep_var(datpath,'data',vals,df=1,tail='',pargs=True,
        mulslice=False,keV=200,NxNy=2**10,repeat=[1,1,1],slice_thick=1,
        opt='srw',fopt='f',ssh='badb')
    return df

def pp_simus():
    df = pp.update_df_info(df_name,files=['patternS'])
    print(df)
    return df

def plot_S():
    df = pd.read_pickle(df_name)
    for obj_pkl in df.index:#[6:]:
        val = df.loc[obj_pkl].data.split('_')[0]
        multi = pp.load_multi_obj(datpath+obj_pkl)
        dutil.plot_npy(multi._outf('patternS'),log=True,
        caxis=[-9,-6],title=r'%s' %val,xylims=0.2,
        name=figpath+val+'log.png',opt='sc')


# gen_simus()
# df = run_simus()
df = pp_simus()
plot_S()
