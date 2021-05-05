#! /bin/bash
#
# Script to merge output ROOT files from all pt-hat bins together

#SUFFIX="2760_PP_Colorless"
#SUFFIX="2760_PbPb_0-5_0.30_2.0_1"
#SUFFIX="2760_PbPb_5-10_0.30_2.0_1"
#SUFFIX="5020_PP_Colorless"
#SUFFIX="5020_PbPb_0-5_0.30_2.0_1"
SUFFIX="5020_PbPb_5-10_0.30_2.0_1"

OUTPUT_DIR_BASE=/rstorage/jetscape/AnalysisResults/499531/v3/$SUFFIX

hadd -f -j 10 $OUTPUT_DIR_BASE/AnalysisResultsFinal.root $OUTPUT_DIR_BASE/Stage0/*/*.root
