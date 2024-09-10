#!/bin/bash

# load environment
source /software/flo/myJETSCAPE/STAT-XSEDE-2021/scripts/venv/bin/activate
pip install uproot
pip install attrs
pip install /software/flo/myJETSCAPE/STAT-XSEDE-2021/scripts
# pip install git+https://github.com/Parsl/parsl.git@refs/pull/3591/merge 

python3 /software/flo/myJETSCAPE/JETSCAPE-analysis/plot/steer_aggregate_and_plot_observables.py