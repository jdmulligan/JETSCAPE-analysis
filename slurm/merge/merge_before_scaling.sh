#! /bin/bash

# Script to merge output ROOT files
#OUTPUT_DIR_BASE="/rstorage/jetscape/AnalysisResults/474425/v2/5020_PP"
OUTPUT_DIR_BASE="/rstorage/jetscape/AnalysisResults/474425/v2/OutputFile_Type5_qhatA10_B100_5020_PbPb_0-10_0.30_2.0_1"

# Loop through pt hat bins and merge files from each pt-hat bin
BINS=(1 2 3 4 5 7 9 11 13 15 17 20 25 30 35 40 45 50 55 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 350 400 450 500 550 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2200 2400 2510)

for ((i=0; i < ${#BINS[@]}-1; ++i))
do
    mkdir -p ${OUTPUT_DIR_BASE}/Stage0/${i}

    FILES=$( find ${OUTPUT_DIR_BASE}/*/${i} -name "*.root" )
    hadd -f -j 10 ${OUTPUT_DIR_BASE}/Stage0/${i}/AnalysisResults.root $FILES
done
