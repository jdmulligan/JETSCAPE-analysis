#!/usr/bin/env python3
import os
import sys
import yaml

sys.path.append('../..')
from jetscape_analysis.analysis import scale_histograms

#suffix = '5020_PP'
suffix = 'OutputFile_Type5_qhatA10_B100_5020_PbPb_0-10_0.30_2.0_1'

base_dir = f'/rstorage/jetscape/AnalysisResults/474425/v2/{suffix}'
pt_hat_dir = f'/rstorage/jetscape/JETSCAPE-AA-events/skim/452210/v2/{suffix}'

config_file = '../../config/TG3.yaml'
with open(config_file, 'r') as stream:
    config = yaml.safe_load(stream)
    pt_hat_bins = config['pt_hat_bins']
    n_pt_hat_bins = len(pt_hat_bins) - 1

for i in range(0, n_pt_hat_bins):
    output_dir = os.path.join(base_dir, f'Stage0/{i}')   

    # Get pt-hat scale factor from file
    pt_hat_min = pt_hat_bins[i]
    pt_hat_max = pt_hat_bins[i+1]
    pt_hat_filename = os.path.join(pt_hat_dir, f'SigmaHardBin{pt_hat_min}_{pt_hat_max}.out')
    with open(pt_hat_filename) as f:
        first_line = f.readline()
        pt_hat_xsec = float(first_line.split(' ')[0])

    scale_histograms.scale_histograms(output_dir, i, pt_hat_xsec, bRemoveOutliers=False)
