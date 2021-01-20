#!/bin/bash
set -e

DIALS_VERSION=$(dials.version | head -1 | awk '{print $NF}')
echo "DIALS version: $DIALS_VERSION"
IFS=. read major minor patch <<< "$DIALS_VERSION"
if [[ "$minor" == "dev" ]]; then
    echo "Development version of DIALS - no version check done"
elif [[ "$major" != "3" || "$minor" -lt "2" ]]; then
    echo "Need DIALS 3.2 or higher"  >&2
    exit 1
fi

# Check script input
if [ "$#" -ne 1 ]; then
    echo "Download IRELOH_ED_Dataset_*.tar.gz from https://zenodo.org/record/1407682#.YAbeFnX7RhE
and unpack into a directory <DATADIR>. Then run ./process <DATADIR>"
    exit 1
fi

PROCDIR=$(pwd)
SCRIPTDIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
DATADIR=$(realpath "$1")
if [ ! -d "$DATADIR" ]; then
    echo "$DATADIR is not found"
    exit 1
fi

# Install/update FormatCBFMiniTimepix
dxtbx.install_format -u\
    https://raw.githubusercontent.com/dials/dxtbx_ED_formats/master/FormatCBFMiniTimepix.py

cat > restraint.phil <<+
refinement
{
  parameterisation
  {
    crystal
    {
      unit_cell
      {
        restraints
        {
          tie_to_target
          {
            values=8.0150,10.015,17.703,90,90,90
            sigmas=0.01,0.01,0.01,0.01,0.01,0.01
          }
        }
      }
    }
  }
}
+

integrate_one () {
    TEMPLATE=$1
    BEAM_CENTRE=$2
    AXIS=$3
    IDX_ASSN_METHOD=$4

    dials.import template="$TEMPLATE"/\
        distance=489\
        slow_fast_beam_centre="$BEAM_CENTRE"\
        geometry.scan.oscillation=0,0.0652\
        goniometer.axis="$AXIS"
    dials.find_spots imported.expt kernel_size=9,9 min_spot_size=8\
        gain=0.8 d_max=15
    dials.index imported.expt strong.refl\
        beam.fix=all goniometer.fix=None detector.fix=distance\
        unit_cell=8.0150,10.015,17.703,90,90,90 space_group=P212121\
        index_assignment.method="$IDX_ASSN_METHOD"
    dials.refine indexed.{expt,refl}\
        scan_varying=False "$PROCDIR"/restraint.phil\
        beam.fix=all goniometer.fix=None detector.fix=distance
    dials.refine refined.{expt,refl} crystal.unit_cell.force_static=true\
        beam.fix=all goniometer.fix=None detector.fix=all
    dials.plot_scan_varying_model refined.expt
    dials.integrate refined.expt refined.refl d_min=0.85
    dials.scale integrated.{expt,refl}
    dials.report scaled.{expt,refl}
}

scale () {
    dials.scale "$PROCDIR"/IRELOH_ED_Dataset_1-dials/scaled.{expt,refl}\
         "$PROCDIR"/IRELOH_ED_Dataset_2-dials/scaled.{expt,refl}\
         "$PROCDIR"/IRELOH_ED_Dataset_3-dials/scaled.{expt,refl}\
         d_min=0.87
    dials.split_experiments scaled.expt scaled.refl
    mv split_0.expt scaled_1.expt
    mv split_1.expt scaled_2.expt
    mv split_2.expt scaled_3.expt
    mv split_0.refl scaled_1.refl
    mv split_1.refl scaled_2.refl
    mv split_2.refl scaled_3.refl
}

# IRELOH_ED_Dataset_1
cd "$PROCDIR"
mkdir -p IRELOH_ED_Dataset_1-dials
cd IRELOH_ED_Dataset_1-dials
TEMPLATE="$DATADIR"/IRELOH_ED_Dataset_1/n14_a004_####.cbf
integrate_one "$TEMPLATE" 266,268 0.998341,-0.0575638,-0.000888014 "local"
cd "$PROCDIR"

# IRELOH_ED_Dataset_2
cd "$PROCDIR"
mkdir -p IRELOH_ED_Dataset_2-dials
cd IRELOH_ED_Dataset_2-dials
TEMPLATE="$DATADIR"/IRELOH_ED_Dataset_2/n14_a006_###.cbf
integrate_one "$TEMPLATE" 269,270 0.99889,-0.0458683,0.0107261 "simple"
cd "$PROCDIR"

# IRELOH_ED_Dataset_3
cd "$PROCDIR"
mkdir -p IRELOH_ED_Dataset_3-dials
cd IRELOH_ED_Dataset_3-dials
TEMPLATE="$DATADIR"/IRELOH_ED_Dataset_3/n14_a009_####.cbf
integrate_one "$TEMPLATE" 269,279 -0.999563,0.0274898,-0.0108268 "simple"
cd "$PROCDIR"

# Scaling
cd "$PROCDIR"
mkdir -p scaling
cd scaling
scale
cd "$PROCDIR"
