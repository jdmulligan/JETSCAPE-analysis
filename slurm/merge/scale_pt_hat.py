#!/usr/bin/env python3
import os
import sys
import yaml

sys.path.append('../..')
from jetscape_analysis.analysis import scale_histograms

#base_dir = '/rstorage/jetscape/AnalysisResults/472497/v2/5020_PP'
base_dir = '/rstorage/jetscape/AnalysisResults/472497/v2/OutputFile_Type5_qhatA10_B100_5020_PbPb_0-10_0.30_2.0_1'

config_file = '../../config/TG3.yaml'
with open(config_file, 'r') as stream:
    config = yaml.safe_load(stream)
    pt_hat_bins = config['pt_hat_bins']
    n_pt_hat_bins = len(pt_hat_bins) - 1

for i in range(0, n_pt_hat_bins):
    output_dir = os.path.join(base_dir, f'Stage0/{i}')   
    scale_histograms.scale_histograms(output_dir, i, bRemoveOutliers=False)
