#!/bin/bash

# load environment
BASEPATH=/software/users/fjonas/myJETSCAPE
source ${BASEPATH}/STAT-XSEDE-2021/scripts/venv/bin/activate
pip install ${BASEPATH}/STAT-XSEDE-2021/scripts
module use /software/users/alice/yasp/software/modules
module load root
python3 ${BASEPATH}/JETSCAPE-analysis/plot/steer_aggregate_and_plot_observables.py
