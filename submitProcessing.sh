#!/bin/bash

# load environment
BASEPATH=/software/users/fjonas/myJETSCAPE
source ${BASEPATH}/STAT-XSEDE-2021/scripts/venv/bin/activate
pip install uproot
pip install attrs
pip install ${BASEPATH}/STAT-XSEDE-2021/scripts
# pip install git+https://github.com/Parsl/parsl.git@refs/pull/3591/merge 
module use /software/users/alice/yasp/software/modules
module load root
python3 ${BASEPATH}/JETSCAPE-analysis/plot/steer_aggregate_and_plot_observables.py
