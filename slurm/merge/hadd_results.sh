#! /bin/bash

# Load modules
module use /software/users/james/heppy/modules
module load heppy/1.0
module list

# Go to dir, and hadd in each dir
cd /rstorage/jetscape/PHYS_RAA/212739/

for d in */ ; do
    echo "$d"
    cd $d
    hadd -j 10 AnalysisResultsFinal.root */*/*.root
    cd ..
done
