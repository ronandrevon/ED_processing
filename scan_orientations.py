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
from libtbx import phil
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util import Sorry, show_mail_handle_errors


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

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    from dials.util.options import OptionParser, flatten_experiments

    usage = "dials.python scan_orientations.py integrated.expt integrated.refl cif_file=mol.cif scan_point=0"

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
    try:
        structure = gemmi.read_small_structure(params.cif_file)
    except Exception as e:
        sys.exit("No structural model provided to rotate")

    # Write out structure at the requested orientation
    orientation = extract_orientation(experiment, params.image)
    basename, ext = os.path.splitext(os.path.basename(params.cif_file))
    filename = basename + ".csv"
    if structure:
        write_rotated(structure, orientation, filename)

    # Find reflections close to the Ewald sphere
    nearby = close_to_Ewald_sphere(reflections, experiments, params.image)


def extract_orientation(exp, image):
    """Extract the crystal orientation at the specified image boundary, i"""
    crystal = exp.crystal
    beam = exp.beam
    scan = exp.scan
    gonio = exp.goniometer

    # Correct for first image number to index scan points
    array_range = scan.get_array_range()
    if image < array_range[0] or image > array_range[1]:
        sys.exit(f"image {i} lies outside the allowed range {array_range}")
    i = image - array_range[0]

    if gonio.num_scan_points > 0:
        S = matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
    else:
        S = matrix.sqr(gonio.get_setting_rotation())

    F = matrix.sqr(gonio.get_fixed_rotation())
    axis = matrix.col(gonio.get_rotation_axis_datum())
    phi = scan.get_angle_from_array_index(i, deg=True)
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

    return SRFU


def write_rotated(structure, orientation, filename):
    sites = structure.get_all_unit_cell_sites()
    tr = gemmi.Transform()
    tr.mat.fromlist(orientation.as_list_of_lists())

    pos = []
    for site in sites:
        pos.append(tr.apply(site.orth(structure.cell)))

    with open(filename, "w") as f:
        f.write("Atom,X,Y,Z,Occ,u_iso\n")
        for s, p in zip(sites, pos):
            f.write(f"{s.element.atomic_number},{p.x},{p.y},{p.z},{s.occ},{s.u_iso}\n")


def close_to_Ewald_sphere(reflections, experiments, scan_point):
    # TODO
    pass


if __name__ == "__main__":
    run()
