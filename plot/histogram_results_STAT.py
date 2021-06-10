"""
  macro for plotting analyzed jetscape events
  """

# This script plots histograms created in the analysis of Jetscape events
#
# Author: James Mulligan (james.mulligan@berkeley.edu)

# General
import os
import sys
import yaml
import argparse
import operator

# Data analysis and plotting
import ROOT
import ctypes
from array import *
import numpy as np
import pandas as pd
import awkward as ak

# Base class
sys.path.append('.')
from jetscape_analysis.base import common_base
from plot import plot_results_STAT_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class HistogramResults(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(HistogramResults, self).__init__(**kwargs)
        self.input_file = input_file
        self.output_dir = output_dir
        
        self.plot_utils = plot_results_STAT_utils.PlotUtils()

        #------------------------------------------------------
        # Read config file
        with open(config_file, 'r') as stream:
            self.config = yaml.safe_load(stream)
            
        self.sqrts = self.config['sqrt_s']
                    
        #------------------------------------------------------
        # Read input file
        self.observables_df = pd.read_parquet(self.input_file)
        #self.observables_df = ak.Array(self.input_file)
        self.weights = self.observables_df['event_weight']
        
        #------------------------------------------------------
        # Read cross-section file
        cross_section_file = 'cross_section'.join(self.input_file.rsplit('observables', 1))
        cross_section_df = pd.read_parquet(cross_section_file)
        self.cross_section = cross_section_df['cross_section'][0]
        self.cross_section_error = cross_section_df['cross_section_error'][0]
        self.n_events = cross_section_df['n_events'][0]
        
        #------------------------------------------------------
        # Create output list to store histograms
        self.output_list = []

        #print(self)
        #print(f'keys: {self.observables_df.keys()}')

    #-------------------------------------------------------------------------------------------
    # Main function
    #-------------------------------------------------------------------------------------------
    def histogram_results(self):
        
        self.histogram_hadron_observables(observable_type='hadron')
        
        self.histogram_hadron_correlation_observables(observable_type='hadron_correlations')
        
        self.histogram_jet_observables(observable_type='inclusive_chjet')
        
        if 'inclusive_jet' in self.config:
            self.histogram_jet_observables(observable_type='inclusive_jet')
            
        if 'semi_inclusive_chjet' in self.config:
            self.histogram_semi_inclusive_chjet_observables(observable_type='semi_inclusive_chjet')
            
        if 'dijet' in self.config:
            self.histogram_jet_observables(observable_type='dijet')
            
        h = ROOT.TH1F('h_n_events', 'h_n_events', 1, 0, 1)
        h.SetBinContent(1, self.n_events)
        self.output_list.append(h)

        self.write_output_objects()

    #-------------------------------------------------------------------------------------------
    # Histogram hadron observables
    #-------------------------------------------------------------------------------------------
    def histogram_hadron_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')
                
        for observable, block in self.config[observable_type].items():
        
            # Construct appropriate binning
            bins = self.plot_utils.bins_from_config(block, observable)
            if not bins.any():
                continue
                
            # Histogram observable
            self.histogram_observable(column_name=f'{observable_type}_{observable}', bins=bins)
            self.histogram_observable(column_name=f'{observable_type}_{observable}_holes', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram hadron correlation observables
    #-------------------------------------------------------------------------------------------
    def histogram_hadron_correlation_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')

        for observable, block in self.config[observable_type].items():

            # Construct appropriate binning
            bins = self.plot_utils.bins_from_config(block, observable)
            if not bins.any():
                continue
                    
            # Histogram observable
            self.histogram_observable(column_name=f'{observable_type}_{observable}', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def histogram_jet_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')

        for observable, block in self.config[observable_type].items():
                
            for jet_R in block['jet_R']:
                print(f'    R = {jet_R}')
                
                # Construct appropriate binning
                bins = self.plot_utils.bins_from_config(block, observable, jet_R=jet_R)
                if not bins.any():
                    continue
                
                if 'SoftDrop' in block:
                    for grooming_setting in block['SoftDrop']:
                        if observable == 'tg_alice' and jet_R == 0.2 and grooming_setting['zcut'] == 0.4:
                            continue
                        else:
                            print(f'      grooming_setting = {grooming_setting}')
                            zcut = grooming_setting['zcut']
                            beta = grooming_setting['beta']
                            column_name = f'{observable_type}_{observable}_R{jet_R}_zcut{zcut}_beta{beta}'
                            self.histogram_observable(column_name=column_name, bins=bins)
                else:
                    self.histogram_observable(column_name=f'{observable_type}_{observable}_R{jet_R}', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram semi-inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def histogram_semi_inclusive_jet_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')
        
        for observable, block in self.config[observable_type].items():

            # Construct appropriate binning
            bins = plot_results_STAT_utils.bins_from_config(block, observable)
            if not bins.any():
                continue
                
            for jet_R in block['jet_R']:
                print(f'    R = {jet_R}')

                column_name = f'{observable_type}_{observable}_R{jetR}_lowTrigger'
                self.histogram_observable(column_name=column_name, bins=bins)
                
                column_name = f'{observable_type}_{observable}_R{jetR}_highTrigger'
                self.histogram_observable(column_name=column_name, bins=bins)
                
                if self.sqrts == '2760':
                    column_name = f'{observable_type}_alice_trigger_pt'
                if self.sqrts == '200':
                    column_name = f'{observable_type}_star_trigger_pt'
                self.histogram_observable(column_name=column_name, bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram a single observable
    #-------------------------------------------------------------------------------------------
    def histogram_observable(self, column_name=None, bins=None):

        # Get column
        col = self.observables_df[column_name]
        
        # Find dimension of observable
        dim_observable = 0
        for i,_ in enumerate(col):
            if len(col[i]) > 0:
                if isinstance(col[i], list):
                    dim_observable = len(col[i][0])
                else:
                    dim_observable = 1
                break
        
        # Construct histogram
        if dim_observable == 1:
            self.histogram_1d_observable(col, column_name=column_name, bins=bins)
        elif dim_observable == 2:
            self.histogram_2d_observable(col, column_name=column_name, bins=bins)
        else:
            return

    #-------------------------------------------------------------------------------------------
    # Histogram a single observable
    #-------------------------------------------------------------------------------------------
    def histogram_1d_observable(self, col, column_name=None, bins=None):

        hname = f'h_{column_name}'
        h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
        h.Sumw2()
        
        # Fill histogram
        for i,_ in enumerate(col):
            if len(col[i]) > 0:
                for value in col[i]:
                    h.Fill(value, self.weights[i])
                    
        # Save histogram to output list
        self.output_list.append(h)
        
    #-------------------------------------------------------------------------------------------
    # Histogram a single observable
    #-------------------------------------------------------------------------------------------
    def histogram_2d_observable(self, col, column_name=None, bins=None):

        hname = f'h_{column_name}'
        h = ROOT.TH2F(hname, hname, len(bins)-1, bins, len(bins)-1, bins)
        h.Sumw2()
        
        # Fill histogram
        for i,_ in enumerate(col):
            if len(col[i]) > 0:
                for value in col[i]:
                    h.Fill(value[0], value[1], self.weights[i])
                    
        # Save histogram to output list
        self.output_list.append(h)

    # ---------------------------------------------------------------
    # Save all ROOT histograms to file
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Save output objects
        output_file = self.input_file.replace('observables', 'histograms').replace('parquet', 'root')
        output_path = os.path.join(self.output_dir, output_file)
        fout = ROOT.TFile(output_path, 'recreate')
        fout.cd()
        for obj in self.output_list:

            types = (ROOT.TH1, ROOT.THnBase)
            if isinstance(obj, types):
                obj.Write()
                obj.SetDirectory(0)
                del obj

        fout.Close()

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    print('Executing histogram_results_STAT.py...')

    # Define arguments
    parser = argparse.ArgumentParser(description='Histogram JETSCAPE observables')
    parser.add_argument(
        '-c',
        '--configFile',
        action='store',
        type=str,
        metavar='configFile',
        default='config/TG3.yaml',
        help='Config file'
    )
    parser.add_argument(
        '-i',
        '--inputFile',
        action='store',
        type=str,
        metavar='inputFile',
        default='observables_5020_0000_00.parquet',
        help='Input file'
    )
    parser.add_argument(
        '-o',
        '--outputDir',
        action='store',
        type=str,
        metavar='outputDir',
        default='/home/jetscape-user/JETSCAPE-analysis/TestOutput',
        help='Output directory for output to be written to'
    )

    # Parse the arguments
    args = parser.parse_args()
    
    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File "{0}" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid inputDir is given, exit
    if not os.path.exists(args.inputFile):
        print('File "{0}" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)
        
    # If output dir doesn't exist, create it
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)

    analysis = HistogramResults(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.histogram_results()
