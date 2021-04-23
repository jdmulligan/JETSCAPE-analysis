#! /bin/bash
#
# Script to merge output ROOT files from all pt-hat bins together
#OUTPUT_DIR_BASE="/rstorage/jetscape/AnalysisResults/474425/v2/5020_PP"
OUTPUT_DIR_BASE="/rstorage/jetscape/AnalysisResults/474425/v2/OutputFile_Type5_qhatA10_B100_5020_PbPb_0-10_0.30_2.0_1"

hadd -f -j 10 $OUTPUT_DIR_BASE/AnalysisResultsFinal.root $OUTPUT_DIR_BASE/Stage0/*/*.root
