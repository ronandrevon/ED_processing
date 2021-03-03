# import importlib as imp
from subprocess import Popen,PIPE
import numpy as np
import binascii,os
from glob import glob
import import_cif
from libtbx import phil
from dials.util.options import OptionParser, reflections_and_experiments_from_files,flatten_experiments,flatten_reflections
from cbflib_adaptbx import compress
from scitbx.array_family import flex
from scitbx import matrix

adxv = os.environ['HOME']+'/bin/adxv.x86_64CentOS7'
Rx180 = np.array([
    [1,0,0 ],
    [0,0,-1],
    [0,1,0 ]])

def gen_xyz(path,name='',cif_file=None,step=1,rep=1,pad=0):
    '''
    - path:postprocess directory containing exp and refl
    - name : prefix to save figures
    - cif_file:full path to cif_file (if None taking first instance of cif_file in directory)
    - step : step for images to take'''
    if not cif_file:cif_file = glob(path+'*.cif')[0]       ;print(cif_file)
    exp,refl    = get_exp_refl(path,'indexed.expt','indexed.refl')
    scan        = exp.scan
    array_range = scan.get_array_range()
    # files = glob(path+'*.cbf')
    # images = [int(f.split('_')[-1].replace('.cbf','')) for f in files];print(images)
    import_cif.import_cif(cif_file,rep=rep)
    for image in range(array_range[0],array_range[1],step):
        i = image - array_range[0]
        orientation = extract_orientation(exp, i)
        filename = path+'%s%s.xyz' %(name,str(i).zfill(4))
        structure = import_cif.import_xyz(cif_file.replace('.cif','.xyz'))
        write_rotated(structure, orientation, filename,pad=pad)

def get_exp_refl(path,exp_file='imported.expt',refl_file='strong.refl'):
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


def write_rotated(structure, orientation, filename,pad=0):
    orientation = orientation.as_mat3()
    pattern,lat_params = structure
    coords = flex.vec3_double(np.array(pattern[:,1:4],dtype=np.double))
    coords = orientation*coords
    # Rotate 180 degrees around x to keep rotation axis
    # coords = np.array(coords)
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
    - im : np.ndarray or str : image to convert to cbf
    - orig_path : the exp data to get the header from
    - out_cbf : name of the ouput cbf
    - pOpt : Show cbf file with adxv if True
    '''
    if isinstance(im,str) :
        im = np.load(im)
        if not out_cbf : out_cbf = in_npy.replace('.npy','.cbf')
    orig_file = glob(orig_path+'*.cbf')[1]; print(orig_file)
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

    im001 = flex.int32(np.array(im*2**15,dtype=np.int32))
    compressed = compress(im001)
    nbytes = len(compressed)

    # Update the header
    pre, post = cbf_header.decode().split(xbsize_record)
    new_xbsize_record = "X-Binary-Size:{0:10d}".format(nbytes)
    if xbsize_record.endswith("\r"):
        new_xbsize_record += "\r"
    new_cbf_header = pre + new_xbsize_record + post

    open(out_cbf, "wb").write(new_cbf_header.encode() + start_tag + compressed + tail)
    # print(colors.green +'file saved : \n' +colors.yellow+out_cbf+colors.black)
    if pOpt:
        p = Popen("%s %s" %(adxv,out_cbf),
            shell=True);p.wait()


# def get_A(exp):
#     crystal_model = exp.crystal
#     cs = crystal_model.symmetry(
#         unit_cell=crystal_model.get_unit_cell(),
#         space_group=crystal_model.get_space_group(),
#     )
#     cb_op = cs.change_of_basis_op_to_reference_setting()
#     cr = crystal_model.change_basis(cb_op)
#     A = matrix.sqr(cr.get_A())
#     return A
#########################################################################################
#### display
#########################################################################################
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
