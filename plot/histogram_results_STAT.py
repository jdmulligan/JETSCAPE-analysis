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

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class HistogramResults(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(HistogramResults, self).__init__(**kwargs)
        self.output_dir = output_dir

        #------------------------------------------------------
        # Read config file
        with open(config_file, 'r') as stream:
            self.config = yaml.safe_load(stream)
            
        self.sqrts = self.config['sqrt_s']
                    
        #------------------------------------------------------
        # Read input file
        self.observables_df = pd.read_parquet(input_file)
        #self.observables_df = ak.Array(input_file)
        self.weights = self.observables_df['event_weight']
        
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

        self.write_output_objects()

    #-------------------------------------------------------------------------------------------
    # Histogram hadron observables
    #-------------------------------------------------------------------------------------------
    def histogram_hadron_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')
                
        for observable, block in self.config[observable_type].items():
        
            # Construct appropriate binning
            bins = self.bins_from_config(block, observable)
            if not bins.any():
                continue
                
            # Histogram observable
            self.histogram_1d_observable(column_name=f'{observable_type}_{observable}', bins=bins)
            self.histogram_1d_observable(column_name=f'{observable_type}_{observable}_holes', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram hadron correlation observables
    #-------------------------------------------------------------------------------------------
    def histogram_hadron_correlation_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')

        for observable, block in self.config[observable_type].items():

            # Construct appropriate binning
            bins = self.bins_from_config(block, observable)
            if not bins.any():
                continue
                    
            # Histogram observable
            self.histogram_1d_observable(column_name=f'{observable_type}_{observable}', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def histogram_jet_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')

        for observable, block in self.config[observable_type].items():

            # Construct appropriate binning
            bins = self.bins_from_config(block, observable)
            if not bins.any():
                continue
                
            for jet_R in block['jet_R']:
                print(f'    R = {jet_R}')
                if 'SoftDrop' in block:
                    for grooming_setting in block['SoftDrop']:
                        if observable == 'tg_alice' and jet_R == 0.2 and grooming_setting['zcut'] == 0.4:
                            continue
                        else:
                            print(f'      grooming_setting = {grooming_setting}')
                            zcut = grooming_setting['zcut']
                            beta = grooming_setting['beta']
                            column_name = f'{observable_type}_{observable}_R{jet_R}_zcut{zcut}_beta{beta}'
                            self.histogram_1d_observable(column_name=column_name, bins=bins)
                else:
                    self.histogram_1d_observable(column_name=f'{observable_type}_{observable}_R{jet_R}', bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram semi-inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def histogram_semi_inclusive_jet_observables(self, observable_type=''):
        print()
        print(f'Histogram {observable_type} observables...')
        
        for observable, block in self.config[observable_type].items():

            # Construct appropriate binning
            bins = self.bins_from_config(block, observable)
            if not bins.any():
                continue
                
            for jet_R in block['jet_R']:
                print(f'    R = {jet_R}')

                column_name = f'{observable_type}_{observable}_R{jetR}_lowTrigger'
                self.histogram_1d_observable(column_name=column_name, bins=bins)
                
                column_name = f'{observable_type}_{observable}_R{jetR}_highTrigger'
                self.histogram_1d_observable(column_name=column_name, bins=bins)
                
                if self.sqrts == '2760':
                    column_name = f'{observable_type}_alice_trigger_pt'
                if self.sqrts == '200':
                    column_name = f'{observable_type}_star_trigger_pt'
                self.histogram_1d_observable(column_name=column_name, bins=bins)

    #-------------------------------------------------------------------------------------------
    # Histogram a single observable
    #-------------------------------------------------------------------------------------------
    def histogram_1d_observable(self, column_name=None, bins=None):

        # Construct histogram
        hname = f'h_{column_name}'
        h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
        h.Sumw2()
        
        # Fill histogram
        col = self.observables_df[column_name]
        for i,_ in enumerate(col):
            if len(col[i]) > 0:
                for value in col[i]:
                    h.Fill(value, self.weights[i])
                    
        # Save histogram to output list
        self.output_list.append(h)

    # ---------------------------------------------------------------
    # Get bin array specified in config block
    # ---------------------------------------------------------------
    def bins_from_config(self, block, observable):
    
        if 'hepdata' in block:
            print(f'  Histogram with hepdata binning for {observable}')
            return self.bins_from_hepdata(block)
        elif 'bins' in block:
            print(f'  Histogram with custom binning for {observable}')
            return np.array(block['bins'])
        else:
            print(f'  Warning: No binning found for {observable}')
            return np.array([])
            
    # ---------------------------------------------------------------
    # Get bin array from hepdata file specified in config block
    # ---------------------------------------------------------------
    def bins_from_hepdata(self, block):

        f = ROOT.TFile(block['hepdata'], 'READ')
        dir = f.Get(block['hepdata_dir'])
        h = dir.Get(block['hepdata_hname'])
        bins = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        return bins

    # ---------------------------------------------------------------
    # Save all ROOT histograms to file
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Fill cross-section
        #self.hCrossSection.SetBinContent(self.pt_hat_bin+1, self.pt_hat_xsec)
        #self.hCrossSection.SetBinError(self.pt_hat_bin+1, self.pt_hat_xsec_err)

        # Set N events
        #self.hNevents.SetBinContent(self.pt_hat_bin+1, self.n_event_max)

        # Save output objects
        outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
        fout = ROOT.TFile(outputfilename, 'recreate')
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
