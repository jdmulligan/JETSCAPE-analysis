#!/bin/bash

# load environment
source /software/flo/myJETSCAPE/STAT-XSEDE-2021/scripts/venv/bin/activate
pip install uproot
pip install /software/flo/myJETSCAPE/STAT-XSEDE-2021/scripts
python3 /software/flo/myJETSCAPE/JETSCAPE-analysis/plot/steer_aggregate_and_plot_observables.py