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
    echo "Download EPICZA_ED_Dataset_*.tar.gz from https://zenodo.org/record/1407682#.YAbeFnX7RhE
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
            values=10.996,12.452,13.218,90,90,90
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
    OSCILLATION=$4
    GAIN=$5

    dials.import template="$TEMPLATE"/\
        distance=489\
        slow_fast_beam_centre="$BEAM_CENTRE"\
        geometry.scan.oscillation="$OSCILLATION"\
        goniometer.axis="$AXIS"
    dials.find_spots imported.expt kernel_size=9,9 min_spot_size=8\
        gain="$GAIN" d_max=15
    dials.index imported.expt strong.refl\
        beam.fix=all goniometer.fix=None detector.fix=distance\
        unit_cell=10.996,12.452,13.218,90,90,90 space_group=P212121
    dials.refine indexed.{expt,refl}\
        scan_varying=False "$PROCDIR"/restraint.phil\
        beam.fix=all goniometer.fix=None detector.fix=distance
    dials.refine refined.{expt,refl}\
        crystal.unit_cell.force_static=true\
        beam.fix=all goniometer.fix=None detector.fix=all
    dials.plot_scan_varying_model refined.expt
    dials.integrate refined.expt refined.refl d_min=0.85
    dials.scale integrated.{expt,refl}
    dials.report scaled.{expt,refl}
}

scale () {
    dials.scale "$PROCDIR"/EPICZA_ED_Dataset_1-dials/scaled.{expt,refl}\
         "$PROCDIR"/EPICZA_ED_Dataset_2-dials/scaled.{expt,refl}\
         "$PROCDIR"/EPICZA_ED_Dataset_3-dials/scaled.{expt,refl}\
         "$PROCDIR"/EPICZA_ED_Dataset_4-dials/scaled.{expt,refl}\
         d_min=0.87
    dials.split_experiments scaled.expt scaled.refl
    mv split_0.expt scaled_1.expt
    mv split_1.expt scaled_2.expt
    mv split_2.expt scaled_3.expt
    mv split_3.expt scaled_4.expt
    mv split_0.refl scaled_1.refl
    mv split_1.refl scaled_2.refl
    mv split_2.refl scaled_3.refl
    mv split_3.refl scaled_4.refl
}

# EPICZA_ED_Dataset_1
cd "$PROCDIR"
mkdir -p EPICZA_ED_Dataset_1-dials
cd EPICZA_ED_Dataset_1-dials
TEMPLATE="$DATADIR"/EPICZA_ED_Dataset_1/n15_a002_####.cbf
integrate_one "$TEMPLATE" 267,269 0.985253,-0.00608263,0.170994 0,0.0662 1
cd "$PROCDIR"

# EPICZA_ED_Dataset_2
cd "$PROCDIR"
mkdir -p EPICZA_ED_Dataset_2-dials
cd EPICZA_ED_Dataset_2-dials
TEMPLATE="$DATADIR"/EPICZA_ED_Dataset_2/n15_a003_####.cbf
integrate_one "$TEMPLATE" 267,268 -0.998318,0.05776,-0.00509128 0,0.0664 1
cd "$PROCDIR"

# EPICZA_ED_Dataset_3
cd "$PROCDIR"
mkdir -p EPICZA_ED_Dataset_3-dials
cd EPICZA_ED_Dataset_3-dials
TEMPLATE="$DATADIR"/EPICZA_ED_Dataset_3/n15_a004_####.cbf
integrate_one "$TEMPLATE" 267,268 0.999264,-0.0181316,-0.0337904 0,0.0664 1
cd "$PROCDIR"

# EPICZA_ED_Dataset_4
cd "$PROCDIR"
mkdir -p EPICZA_ED_Dataset_4-dials
cd EPICZA_ED_Dataset_4-dials
TEMPLATE="$DATADIR"/EPICZA_ED_Dataset_4/n15_a005_####.cbf
integrate_one "$TEMPLATE" 267,268 0.999821,-0.0187882,-0.00243261 0,0.071 1
cd "$PROCDIR"

# Scaling
cd "$PROCDIR"
mkdir -p scaling
cd scaling
scale
cd "$PROCDIR"
