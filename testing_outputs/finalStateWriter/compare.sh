#!/usr/bin/env bash

comparisonDir=testing_outputs/finalStateWriter

cd ../../

source .venv/bin/activate

python3 jetscape_analysis/analysis/reader/skim_ascii.py -i /home/rehlers/code/jetscape/stat-xsede-2021/containers/test_build/jetscape_stat/build/test_out_final_state_hadrons.dat -o $comparisonDir/test_out_final_state_hadrons.parquet
python3 jetscape_analysis/analysis/reader/skim_ascii.py -i /home/rehlers/code/jetscape/stat-xsede-2021/containers/test_build/jetscape_stat/build/test_out_final_state_partons.dat -o $comparisonDir/test_out_final_state_partons.parquet
python3 jetscape_analysis/analysis/reader/skim_ascii.py -i /home/rehlers/code/jetscape/stat-xsede-2021/containers/test_build/jetscape_stat/build/test_out_final_state_hadrons_reference.dat -o $comparisonDir/test_out_final_state_hadrons_reference.parquet
python3 jetscape_analysis/analysis/reader/skim_ascii.py -i /home/rehlers/code/jetscape/stat-xsede-2021/containers/test_build/jetscape_stat/build/test_out_final_state_partons_reference.dat -o $comparisonDir/test_out_final_state_partons_reference.parquet

cd -

python3 compare.py

