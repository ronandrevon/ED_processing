import importlib as imp
import numpy as np, matplotlib.pyplot as plt
import binascii,os,optparse
from subprocess import Popen,PIPE
from glob import glob
from libtbx import phil
from dials.util.options import OptionParser, reflections_and_experiments_from_files,flatten_experiments,flatten_reflections
from cbflib_adaptbx import compress
from scitbx.array_family import flex
from scitbx import matrix
import import_cif; imp.reload(import_cif)
from utils import glob_colors as colors

adxv = os.environ['HOME']+'/bin/adxv.x86_64CentOS7'


def gen_xyz(exp_path,name='',cif_file=None,step=1,rep=1,pad=0,image=-1,i=-1,bfacts=None,popt=0):
    ''' Generate a set of .xyz from cif_file with the same orientations as an experiment
    - exp_path  : Experiment directory containing the .exp and .refl files
    - name      : prefix to save the files as <prefix><num_image>.xyz
    - cif_file  : full path to cif_file (if None taking first instance of cif_file in directory)
    - rep       : list3 or int : crystal size will be rep[0] x rep[1] x rep[2]
    - pad       : amount of padding around the crystal in super cell units
    - image/i   : int>=0 - Still image mode : image number/index number
    - step      : images[0:-1:step] will be used to generate the .xyz files
    - popt      : show structures in xy plane if True
    Returns :
    - saves the .xyz files
    '''
    if not cif_file:cif_file = glob(exp_path+'*.cif')[0]              #;print(cif_file)
    experiments,refl = get_exp_refl(exp_path,'indexed.expt','indexed.refl')
    exp         = experiments[0]
    scan        = exp.scan
    array_range = scan.get_array_range()
    import_cif.import_cif(cif_file,rep=rep,bfacts=bfacts)
    def rotate_structure(image):
        i = image - array_range[0]
        orientation = extract_orientation(exp, i)
        filename = '%s%s.xyz' %(name,str(image).zfill(4))
        structure = import_cif.import_xyz(cif_file.replace('.cif','.xyz'))
        write_rotated(structure, orientation, filename,pad=pad)
        if popt:fig,ax = show_xyz(filename,opts='xy');plt.show()

    if i>=0:
        rotate_structure(i+array_range[0])
    elif image>=0:
        rotate_structure(image)
    else:
        for image in np.arange(array_range[0],array_range[1],step):
            rotate_structure(image)


def get_exp_refl(path,exp_file='imported.expt',refl_file='strong.refl'):
    '''Get experiments and reflections[0] form an .expt and .refl files
    Returns
    - experiments,reflections[0]
    '''
    if not refl_file:refl_file=exp_file.replace('.expt','.refl')
    phil_scope  = phil.parse(""" """)
    parser      = OptionParser(usage='',phil=None,read_experiments=True,read_reflections=True)
    args        = [path+exp_file,path+refl_file]
    params, _ = parser.parse_args(args, show_diff_phil=True)
    reflections, experiments = reflections_and_experiments_from_files(params.input.reflections,params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)
    exp  = experiments[0]
    refl = reflections[0]
    return experiments,refl

def extract_orientation(exp, i):
    """Extract the crystal orientation at the specified image boundary, i"""
    crystal = exp.crystal
    beam    = exp.beam
    scan    = exp.scan
    gonio   = exp.goniometer              #; print(gonio.get_rotation_axis())


    if gonio.num_scan_points > 0:
        S = matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
    else:
        S = matrix.sqr(gonio.get_setting_rotation())

    F = matrix.sqr(gonio.get_fixed_rotation())
    axis = matrix.col(gonio.get_rotation_axis_datum())
    phi = scan.get_angle_from_array_index(i, deg=True)      #;print(phi)
    R = matrix.sqr(axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True))

    if crystal.num_scan_points > 0:
        UB = matrix.sqr(crystal.get_A_at_scan_point(i))
        U = matrix.sqr(crystal.get_U_at_scan_point(i))
    else:
        UB = matrix.sqr(crystal.get_A())
        U = matrix.sqr(crystal.get_U())

    # Construct full setting matrix for the scan point i
    # SRFUB = S * R * F * UB
    SRFU = S * R * F * U
    return SRFU


Rx180 = np.array([
    [1,0,0 ],
    [0,0,-1],
    [0,1,0 ]])
def write_rotated(structure, orientation, filename,pad=0):
    orientation = orientation.as_mat3()
    pattern,lat_params = structure
    coords = flex.vec3_double(np.array(pattern[:,1:4],dtype=np.double))
    coords = orientation*coords
    # coords = np.array(coords)
    #### Rotate 180 degrees around x to keep rotation axis
    coords = Rx180.dot(np.array(coords).T).T
    coords, lat_params = apply_padding(coords,lat_params,pad)
    pattern[:,1:4] = coords
    import_cif.make_xyz(filename,pattern,lat_params,fmt='%.4f')


def apply_padding(coords,lat_params,pad):
    ax = coords[:,0].max()-coords[:,0].min()
    by = coords[:,1].max()-coords[:,1].min()
    cz = coords[:,2].max()-coords[:,2].min()
    coords[:,0] += ax*pad-coords[:,0].min()
    coords[:,1] += by*pad-coords[:,1].min()
    coords[:,2] -= coords[:,2].min()
    lat_params[0] = ax*(1+2*pad)
    lat_params[1] = by*(1+2*pad)
    lat_params[2] = cz #*(1+2*pad)
    return coords, lat_params


def save_cbf(im,orig_path,out_cbf=None,pOpt=False):
    ''' convert np.ndarray to cbf
    - im        : np.ndarray or str : image to convert to cbf
    - orig_path : the path to exp data to get the header from those images
    - out_cbf   : name of the .cbf output file
    - pOpt      : Show cbf file with adxv if True (remember to set the path to adxv l.13)
    '''
    from utils import glob_colors as colors
    if isinstance(im,str) :
        im = np.load(im)
        if not out_cbf : out_cbf = in_npy.replace('.npy','.cbf')
    orig_file = glob(orig_path+'*.cbf')[1] #;print(orig_file)
    # orig_file = "/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/n14_a004_0484.cbf"
    start_tag = binascii.unhexlify("0c1a04d5")

    with open(orig_file,'rb') as cbf : data = cbf.read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[: data_offset - 4]

    fast   = 0
    slow   = 0
    length = 0
    for record in cbf_header.decode().split("\n"):
        if "X-Binary-Size-Fastest-Dimension" in record:
            fast = int(record.split()[-1])
        elif "X-Binary-Size-Second-Dimension" in record:
            slow = int(record.split()[-1])
        elif "X-Binary-Size:" in record:
            xbsize_record = record
            length = int(record.split()[-1])

    tail = data[data_offset + length :]

    # print(im.shape,fast,slow,length)
    im001 = flex.int32(im) #np.array(im,dtype=np.int32))
    compressed = compress(im001)
    nbytes = len(compressed)

    # Update the header
    pre, post = cbf_header.decode().split(xbsize_record)
    new_xbsize_record = "X-Binary-Size:{0:10d}".format(nbytes)
    if xbsize_record.endswith("\r"):
        new_xbsize_record += "\r"
    new_cbf_header = pre + new_xbsize_record + post

    open(out_cbf, "wb").write(new_cbf_header.encode() + start_tag + compressed + tail)
    print(colors.green +'file saved : \n' +colors.yellow+out_cbf+colors.black)
    if pOpt:
        p = Popen("%s %s" %(adxv,out_cbf),
            shell=True);p.wait()

def convert_I(im,cap=2**16,m=2**30):
    ''' converts a float intensity image to uint32
    - im  : str - image as a the .npy file
    - m   : int - multiplying factor
    - cap : int - max value to apply a cap
    '''
    qx,qy,I = np.load(im)
    I/=I.max()
    I = np.array(I*m,dtype=np.uint32)
    I[I>cap] = cap
    # print(I.min(),I.max())
    return qx,qy,I

def npy2cbf(path,image,expdata_path,cap=2**15,m=2**32):
    npy_file  = get_image(path,image,fmt='_autoslic_patternS.npy')
    base_file = os.path.basename(npy_file).split('_')[0]
    cbf_file  = path+'%s_001_%s.cbf' %(base_file,str(image).zfill(4))
    qx,qy,I   = convert_I(npy_file,cap=cap,m=m)
    save_cbf(I,expdata_path,cbf_file,pOpt=0)
    # return qx,qy,I


#########################################################################################
#### display
#########################################################################################
def plot_npy(npy_file,title=None,log=False,**kwargs):
    from utils import displayStandards as dsp   #;imp.reload(dsp)
    qx,qy,I = np.load(npy_file)
    # print('npy : Imax=%d, Imean=%d' %(I.max(),I.mean()))
    I/=I.max()
    if log:
        I[I<1e-10]=1e-10
        I = np.log10(I)
    if not title:title='%s' %os.path.basename(npy_file).replace('_',' '),
    dsp.stddisp(im=[qx,-qy,I],
        title=title,
        imOpt='cv',axPos='V',cmap='binary',**kwargs)

def plot_cbf(cbf_file,caxis=[0,100],**kwargs):
    import multislice.mupy_utils as mut         #;imp.reload(mut)
    imV = mut.Image_viewer(cbf_file,sym=0,pad=4)
    print('cbf : Imax=%d, Imean=%d' %(imV.image.max(),imV.image.mean()))
    return imV.show_image(stack=0,lab='p',cmap='binary',
        caxis=caxis,imOpt='c', axPos='V',**kwargs)

def show_xyz(xyz_name,**kwargs):
    import multislice.mupy_utils as mut #; imp.reload(mut)
    return mut.show_grid(xyz_name,**kwargs)

def xyz_gif(dpath,rpath,opath,name,images=None,xylims=None):
    '''
    dpath : data path
    rpath : all figures path
    opath : path to save the .gif
    images : images to save to fig
    '''
    import multislice.mupy_utils as mut ; imp.reload(mut)
    cmd = "cd %s;for d in $(ls);do rm $d/*;done" %(rpath)
    p = Popen(cmd, shell=True,stderr=PIPE,stdout=PIPE);p.wait()
    if not images:images = np.sort([int(f.split('_')[-1].replace('.xyz','')) for f in glob(dpath+name+'*.xyz')])
    for i in images :
        print(i)
        i_str = str(i).zfill(4)
        xyz_name = dpath+'%s%s.xyz' %(name,i_str)
        mut.show_grid(xyz_name,'xy',name=rpath+'xy/%s.png' %(i_str),opt='sc',xylims=xylims['xy'])
        mut.show_grid(xyz_name,'yz',name=rpath+'yz/%s.png' %(i_str),opt='sc',xylims=xylims['yz'])
        mut.show_grid(xyz_name,'xz',name=rpath+'xz/%s.png' %(i_str),opt='sc',xylims=xylims['xz'])

    cmd = '''cd %s;
    for d in $(ls);do
        gif_file=%s$d.gif;echo $gif_file
        convert -delay 20 -loop 0 $d/*.png $d/$gif_file
        cp $d/$gif_file %s/$gif_file
    done
    ''' %(rpath,name,opath)
    # print(cmd)
    p = Popen(cmd, shell=True,stderr=PIPE,stdout=PIPE)
    p.wait()
    o,e = p.communicate();
    print(o,e)
    # return p

# def show_phi(path):
#     from utils import displayStandards as dsp
#     exp,refl = get_exp_refl(path,'indexed.expt','indexed.refl')
#     scan    = exp.scan
#     i1,i2   = scan.get_array_range()
#     n       = i2-i1
#     phi     = [scan.get_angle_from_array_index(i, deg=True) for i in range(n)]
#     dsp.stddisp([range(n),phi,'b'],labs=['i','$\phi(deg)$'])


####################################################################################
##### post process parser
####################################################################################
def pp_parser():
    usage = '''\n
python3 ireloh_pp.py [--dir=<data_dir>] [--image=<image>] [--opts=<str>] [--plot_opts=<str>] [-h]\n
example :\n
run ireloh_pp.py -d'simus/ireloh2/ireloh2_pad2/' -e'simus/ireloh2/exp/' -i0 o's' -p'pP'
    '''
    description = '''Post processing simulated data'''
    parser = optparse.OptionParser(usage=usage,description=description)#,epilog=epilog,)

    base_dir = 'simus/dat/ireloh/ireloh2/'
    exp_dir = base_dir+'exp/'
    sim_dir = base_dir+'ireloh2_pad2/'
    parser.add_option("-d","--dir",default=sim_dir,action="store",type="string",dest="sim_dir",
        help="str - full path to simulation data")
    parser.add_option("-e","--exp",default=exp_dir,action="store",type="string",dest="exp_dir",
        help="str - full path to experimental data(for cbf header and comparison purposes )")
    parser.add_option("-i","--image",default=-1,action="store",type="int",dest="image",
        help="int - still image number to work with (it will try the image index if it cannot find the actual number)")
    parser.add_option("-o","--opts",default="",action="store",type="string",dest="opts",
        help="str - s(save to cbf) i(dials.import) f(dials.find_spots) h(miller index)")
    parser.add_option("-p","--plot_opts",default="",action="store",type="string",dest="plot_opts",
        help="str - p(plot the .npy) P(plott the .cbf) A(use adxv) I(use dials.image_viewer)")
    return parser

def get_image(path,image,fmt='.cbf'):
    cbf_files = glob(path+'*%s' %fmt)#;print(cbf_files)
    cbf_base  = '_'.join(cbf_files[0].replace(fmt,'').split('_')[:-1])  #;print(cbf_base)
    str_image = str(image).zfill(4)
    cbf_file  = '%s_%s%s' %(cbf_base,str_image,fmt)       #;print(cbf_file)
    if not os.path.exists(cbf_file):
        if len(cbf_files)<=image:
            raise Exception('image %s not found, Available images :\n%s' %(cbf_file,'\n'.join(cbf_files)) )
        warn = '''warning:
        image not found :\n%s
        using instead image :\n%s
        ''' %(colors.yellow+cbf_file+colors.black,colors.green+cbf_files[image]+colors.black )
        print(warn)
        cbf_file  = cbf_files[image]
    return cbf_file
