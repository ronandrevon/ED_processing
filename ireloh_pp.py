'''post processing simulated data for ireloh'''
import os,sys
import numpy as np
from glob import glob
from subprocess import Popen
from libtbx import phil
from scitbx.array_family import flex
from scitbx import matrix
from dials.command_line import find_spots
import dutil
import importlib as imp;imp.reload(dutil)
# interactive_console(); 1/0 #XXXXX DEBUG

parser = dutil.pp_parser()
options,args = parser.parse_args()
path,expdata_path,image,opts,plot_opts = options.__dict__.values()#;print(options,args)
path+='/';expdata_path+='/'



still_image=image>=0
cbf_file=''
if still_image:
    cbf_file=dutil.get_image(path,image,fmt='.cbf')
    npy_file=dutil.get_image(path,image,fmt='_autoslic_patternS.npy')

# def for_all(path,func):
#     files = [os.path.basename(f) for f in glob(path+'*.xyz')]
#     for i,f in enumerate(files):
#         image = int(f.split('_')[-1].replace('.xyz',''))
#         func(image)


if 's' in opts:
    if still_image:dutil.npy2cbf(path,image,expdata_path,cap=2**15,m=2**32)
    else:
        files = [os.path.basename(f) for f in glob(path+'*_autoslic_patternS.npy')]
        for i,f in enumerate(files):
            image = int(f.replace('_autoslic_patternS.npy','').split('_')[-1])
            dutil.npy2cbf(path,image,expdata_path,cap=2**15,m=2**32)

if 'i' in opts:
    cbfs = ['*.cbf',os.path.basename(cbf_file)][still_image]
    cmd = '''
    cd %s
    dials.import %s
    '''%(path,cbfs)#;print(cmd)
    p = Popen(cmd, shell=True); p.wait()

if 'f' in opts:
    cbfs = ['imported.expt',os.path.basename(cbf_file)][still_image]
    p = Popen('''
    cd %s
    dials.find_spots nproc=4 kernel_size=9,9 min_spot_size=8 gain=0.8 d_max=15 threshold=60 %s '''
    %(path,cbfs), shell=True); p.wait()
    
    p = Popen('''
        cd %s &&
        dials.image_viewer imported.expt strong.refl brightness=50 \
        show_shoebox=True show_ctr_mass=True show_miller_indices=True \
        show_predictions=False show_max_pix=False show_threshold_pix=False show_all_pix=False
        ''' %path,
        shell=True); p.wait()

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



#################################################################################
#### display
#################################################################################
if 'I' in plot_opts and 'f' not in opts:
    im = ['imported.expt strong.refl',os.path.basename(cbf_file)][still_image]
    p = Popen('''
        cd %s &&
        dials.image_viewer %s brightness=40 \
        show_shoebox=True show_ctr_mass=True show_miller_indices=True \
        show_predictions=False show_max_pix=False show_threshold_pix=False show_all_pix=False
        ''' %(path,im),
        shell=True); p.wait()

if 'p' in plot_opts and still_image:
    dutil.plot_npy(npy_file,caxis=[0,1e-5],xylims=0.3)

if 'P' in plot_opts and still_image:
    dutil.plot_cbf(cbf_file,caxis=[0,100])#,xylims=)
