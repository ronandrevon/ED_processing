import importlib as imp
import os
import numpy as np
from subprocess import Popen
from libtbx import phil
from scitbx.array_family import flex
from scitbx import matrix
import dutil                                ;imp.reload(dutil)
from dials.command_line import find_spots
# interactive_console(); 1/0 #XXXXX DEBUG


expdata_path = '/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_2/'
simdata_path = '/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/ireloh/ireloh2/'

opts ='' #s(save to cbf) i(dials.import) f(dials.find_spots) h(miller index)
plot_opts='P' #p(plot_npy) P(plt.image_viewer) A(adxv) I(idials.image_viewer)


# work with simulated data
image = 600 # the image number to work with
path  = simdata_path
str_image = str(image).zfill(4)
im        = 'ireloh2_%s_autoslic_patternS.npy' %str_image
cbf_file  = simdata_path+'ireloh2_001_%s.cbf' %str_image


# work with exp data
path = os.path.dirname(__file__)+'/IRELOH_ED_Dataset_2-dials/' #;print(path)
cbf_file  = expdata_path+'n14_a006_160.cbf'


## save file
def convert_I(im):
    qx,qy,I  = np.load(simdata_path+im)
    I/=I.max()
    I = np.array(I*2**30,dtype=np.uint32)
    cap = 2**16
    I[I>cap] = cap
    return qx,qy,I

if 's' in opts:
    qx,qy,I = convert_I(im)
    dutil.save_cbf(I,expdata_path,cbf_file,pOpt=0)

if 'i' in opts:
    p = Popen('''
    cd %s
    dials.import %s'''
    %(cbf_file), shell=True); p.wait()

if 'f' in opts:
    p = Popen('''
    cd %s
    dials.find_spots nproc=4 kernel_size=9,9 min_spot_size=8 gain=0.8 d_max=15 threshold=60 %s '''
    %(cbf_file), shell=True); p.wait()

if 'h' in opts:
    experiments,refl = dutil.get_exp_refl(path,'imported.expt','strong.refl')
    exp = experiments[0]
    imageset = exp.imageset
    detector = imageset.get_detector()
    det = detector[0]

    refl.centroid_px_to_mm(experiments)
    refl.map_centroids_to_reciprocal_space(experiments,calculated=False,crystal_frame=True)
    xy_px = np.array([np.array(r)   for r in refl["xyzobs.px.value"]  ])#.select(sel).parts() ])
    xy_mm = np.array([np.array(r)   for r in refl["xyzobs.mm.value"]  ])#.select(sel).parts() ])
    xy    = [flex.vec2_double(np.array([x[:2]],dtype=np.double)) for x in xy_mm]
    xy_s1 = np.array([det.get_lab_coord(x[0]) for x in xy  ])
    rlp   = np.array([np.array(r)   for r in refl["rlp"]  ])#.select(sel).parts()])
    print('xy_px : ',xy_px)
    print('xy_mm : ',xy_mm)
    print('xy_s1 : ',xy_s1)
    print('rlp   : ',rlp  )

    experiments,refl = dutil.get_exp_refl(path,'indexed.expt','indexed.refl')
    exp = experiments[0]
    # A = dutil.get_A(exp)
    A = exp.crystal.get_A()
    A = np.array(A).reshape((3,3))
    print('checking indexing formula')
    Ainv = np.linalg.inv(A)

    hkl0 = Ainv.dot(rlp.T).T
    #from miller indices
    hkl      = np.array([np.array(r)   for r in refl["miller_index"]])
    print('miller index   :',hkl)
    print('hkl from Ainv  :',hkl0 )


if 'I' in plot_opts:
    p = Popen('''
        cd %s &&
        dials.image_viewer imported.expt strong.refl brightness=100 \
        show_shoebox=True show_ctr_mass=True show_miller_indices=True \
        show_predictions=False show_max_pix=False show_threshold_pix=False show_all_pix=False
        ''' %simdata_path,
        shell=True); p.wait()


def plot_npy(im):
    from utils import displayStandards as dsp   #;imp.reload(dsp)
    rpath = '/home/tarik/Documents/git/ccp4/reports/ireloh/'
    qx,qy,I  = convert_I(im) #;print(I.shape)
    print('Before cbf : Imax=%d, Imean=%d' %(I.max(),I.mean()))
    dsp.stddisp(im=[qx,-qy,I],
        caxis=[0,28000],title=im,#,rings=[0.5,1,2],lw=2,
        xylims=0.3,imOpt='cv',axPos='V',cmap='binary',
        opt='ps',name=rpath+'before.png')
if 'p' in plot_opts:
    plot_npy(im)

if 'P' in plot_opts:
    import multislice.mupy_utils as mut         #;imp.reload(mut)
    rpath = '/home/tarik/Documents/git/ccp4/reports/ireloh/'
    imV = mut.Image_viewer(cbf_file,sym=0)
    print('After  cbf : Imax=%d, Imean=%d' %(imV.image.max(),imV.image.mean()))
    fig,ax = imV.show_image(stack=0,caxis=[0,100],lab='p',cmap='binary',
        imOpt='c', axPos='V',#xylims=[10,18,10,18],
        opt='p',name=rpath+'after.png')

    # ax.scatter(xy_mm[:,0],xy_mm[:,1],50,['r']*xy_mm.shape[0], marker='+')
    # fig.show()
    # scat= [xy_mm[0],xy_mm[1],'b']
    # dsp.stddisp(fig=fig,ax=ax,scat=scat,labs=['$x(mm)$','$y(mm)$'],ms=50,caxis=[0,5000])
