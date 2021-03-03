"""
Rotate a structural model to the orientation given by a DIALS experiment at
a particular scan-point. Output the rotated model plus a list of reflections
close to the Ewald sphere at that orientation, with their extinction errors

Example::

  dials.python scan_orientations.py integrated.expt integrated.refl\
      cif_file=mol.cif image=0
"""

from __future__ import absolute_import, division, print_function

import sys, os
import gemmi
from scitbx import matrix
from scitbx.array_family import flex
from libtbx import phil
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util import Sorry, show_mail_handle_errors
import pytest
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import import_cif

phil_scope = phil.parse(
    """\
cif_file = None
    .type = path
    .help = "File containing the small molecule CIF structural model to rotate"

image = 0
    .type = int
    .help = "Image boundary number for which to extract the crystal orientation."
            "This is a 0-based index, such that if the first image is"
            "image_0001, then 0 here refers to orientation at the start"
            "that image and 1 refers to the orientation at the end of"
            "that image, which is the boundary between image_0001 and"
            "image_0002."

test = False
    .type = bool
    .help = "Run in-line tests during operation"
""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    from dials.util.options import OptionParser, flatten_experiments

    usage = "dials.python scan_orientations.py integrated.expt integrated.refl cif_file=mol.cif image=0"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        epilog=__doc__,
    )

    params, _ = parser.parse_args(args, show_diff_phil=True)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Try to load the integrated geometry models and reflection data
    nexp = len(experiments)
    if nexp == 0 or len(reflections) == 0:
        parser.print_help()
        return
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]
    if len(experiments) > 1:
        sys.exit("Only one experiment can be imported at present")
    experiment = experiments[0]

    # Try to load the structural model
    # try:
        # structure = gemmi.read_small_structure(params.cif_file)
    # except Exception as e:
    #     sys.exit("No structural model provided to rotate")
    structure = import_cif.import_xyz(params.cif_file)

    # Write out structure at the requested orientation
    orientation = extract_orientation(experiment, params.image)
    basename, ext = os.path.splitext(os.path.basename(params.cif_file))
    filename = params.cif_file.replace('.xyz','')+'_rotated%d.xyz' %params.image  #basename + ".csv"
    # if structure:
    if params.test:
        test_rotation(structure, orientation)
    new_pos = write_rotated(structure, orientation, filename,pad=0.5)

    #if params.test:
    #    plot_atom_positions(new_pos)

    # Find reflections close to the Ewald sphere
    nearby = close_to_Ewald_sphere(reflections, experiments, params.image)


def extract_orientation(exp, image):
    """Extract the crystal orientation at the specified image boundary, i"""
    crystal = exp.crystal
    beam = exp.beam
    scan = exp.scan
    gonio = exp.goniometer
    print(gonio.get_rotation_axis())
    # Correct for first image number to index scan points
    array_range = scan.get_array_range()
    if image < array_range[0] or image > array_range[1]:
        sys.exit(f"image {image} lies outside the allowed range {array_range}")
    i = image - array_range[0]

    if gonio.num_scan_points > 0:
        S = matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
    else:
        S = matrix.sqr(gonio.get_setting_rotation())

    F = matrix.sqr(gonio.get_fixed_rotation())
    axis = matrix.col(gonio.get_rotation_axis_datum())
    phi = scan.get_angle_from_array_index(i, deg=True);print(phi)
    R = matrix.sqr(axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True))

    if crystal.num_scan_points > 0:
        UB = matrix.sqr(crystal.get_A_at_scan_point(i))
        U = matrix.sqr(crystal.get_U_at_scan_point(i))
    else:
        UB = matrix.sqr(crystal.get_A())
        U = matrix.sqr(crystal.get_U())

    # Construct full setting matrix for the scan point i
    SRFUB = S * R * F * UB
    SRFU = S * R * F * U

    # SFRUB is the orthogonalisation matrix for the reciprocal space laboratory
    # frame. We want the real space fractionalisation matrix, which is its
    # transpose (https://dials.github.io/documentation/conventions.html)
    frac_mat = SRFUB.transpose()

    # Now get the real space orthogonalisation matrix to calculate the real
    # space cell vectors
    orthog_mat = frac_mat.inverse()
    h = matrix.col((1, 0, 0))
    k = matrix.col((0, 1, 0))
    l = matrix.col((0, 0, 1))
    a = orthog_mat * h
    b = orthog_mat * k
    c = orthog_mat * l
    print(f"At image boundary {image}, the unit cell vectors are")
    print(f"a: {a[0]:.6f}, {a[1]:.6f} {a[2]:.6f}")
    print(f"b: {b[0]:.6f}, {b[1]:.6f} {b[2]:.6f}")
    print(f"c: {c[0]:.6f}, {c[1]:.6f} {c[2]:.6f}")


    # SRFU =  U
    return SRFU


def write_rotated(structure, orientation, filename,pad=0):
    # sites = structure.get_all_unit_cell_sites()
    # SFRU=matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(45,deg=True)
    orientation = orientation.as_mat3()

    # tr = gemmi.Transform()
    # tr.mat.fromlist(orientation.as_list_of_lists())
    pattern,lat_params = structure
    coords = flex.vec3_double(np.array(pattern[:,1:4],dtype=np.double))
    coords = orientation*coords
    # TODO
    # Rotate 180 degrees around x to keep rotation axis

    # Rx = np.array([
    #     [1,0,0 ],
    #     [0,0,-1],
    #     [0,1,0 ]])
    # coords = Rx.dot(np.array(coords).T).T
    # coords, lat_params = apply_padding(coords,lat_params,pad)
    pattern[:,1:4] = coords
    import_cif.make_xyz(filename,pattern,lat_params,fmt='%.4f')

    # pos = []
    # sites = coords
    # for site in sites:
    #     site = flex.vec3_double(np.array(site,dtype=np.double))
    #     # pos.append(tr.apply(site.orth(structure.cell)))
    #     pos.append(orientation*coords)
    # with open(filename, "w") as f:
        # f.write("Atom,X,Y,Z,Occ,u_iso\n")
        # for s, p in zip(sites, pos):
        #     f.write(f"{s.element.atomic_number},{p.x},{p.y},{p.z},{s.occ},{s.u_iso}\n")

    # return pos

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

def plot_atom_positions(positions):
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    X = []
    Y = []
    Z = []
    for pos in positions:
        X.append(pos.x)
        Y.append(pos.y)
        Z.append(pos.z)

    ax.scatter(X, Y, Z)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    #ax.set_aspect("equal")
    set_axes_equal(ax)

    plt.show()
    return

def test_rotation(structure, orientation):
    """Test gemmi rotation transformation does what we expect"""
    sites = structure.get_all_unit_cell_sites()
    tr = gemmi.Transform()
    tr.mat.fromlist(orientation.as_list_of_lists())

    orig_pos = [site.orth(structure.cell) for site in sites]
    plot_atom_positions(orig_pos)
    new_pos = [tr.apply(pos) for pos in orig_pos]

    for vector, rotated in zip(orig_pos, new_pos):
        rotated2 = orientation * matrix.col(vector.tolist())
        for i in range(3):
            assert rotated2[i] == pytest.approx(rotated[i])


def close_to_Ewald_sphere(reflections, experiments, scan_point):
    # TODO
    pass


if __name__ == "__main__":
    run()
    # import_cif.show_grid('ireloh_rotated484.xyz',opts='xy',ms=1,opt='')
    # plt.show()
