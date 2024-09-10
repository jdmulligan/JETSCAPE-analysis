#! /bin/bash

# This script should be run to initialize the analysis environment
# if you are using the JETSCAPE docker container.

# Load heppy module
source /usr/local/init/profile.sh
source /software/flo/myJETSCAPE/STAT-XSEDE-2021/scripts/venv/bin/activate
module load heppy/1.0
