#!/usr/bin/env python3

""" Base class to analyze a JETSCAPE output file: do jet-finding and produce a ROOT file for each pt-hat bin.

You should create a user class that inherits from this one. See analyze_events_PHYS.py for an example.

The outputdir should contain the JETSCAPE output files in the structure generated by the PHYS WG output structure.

See README for pre-requisites.

.. codeauthor:: James Mulligan <james.mulligan@berkeley.edu>, UC Berkeley
"""

from __future__ import print_function

# General
import os
import yaml

# Analysis
import itertools
import ROOT
import pandas as pd
import numpy as np

# Fastjet via python (from external library heppy)
import fjext

from jetscape_analysis.analysis import scale_histograms
from jetscape_analysis.base import common_base

################################################################
class AnalyzeJetscapeEvents_BasePHYS(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file="", input_file="", output_dir="", **kwargs):
        super(AnalyzeJetscapeEvents_BasePHYS, self).__init__(**kwargs)
        self.config_file = config_file

        self.input_file_hadrons = input_file
        self.input_file_partons = ''

        # Get pt-hat from filename
        filename = self.input_file_hadrons.split('/')[-1]
        suffix = filename.split('Bin')[1]
        self.pt_hat_min = int(suffix.split('_')[0])
        self.pt_hat_max = int(suffix.split('_')[1])
        self.index = suffix.split('_')[2]
        self.output_dir = os.path.join(output_dir, self.index)

        # Get pt-hat scale factor from file in same directory
        self.input_dir = os.path.dirname(input_file)
        pt_hat_filename = os.path.join(self.input_dir, '../SigmaHardBin{}_{}.out'.format(self.pt_hat_min, self.pt_hat_max))
        with open(pt_hat_filename) as f:
            first_line = f.readline()
            self.pt_hat_xsec = float(first_line.split(' ')[0])
            self.pt_hat_xsec_err = float(first_line.split(' ')[1])

        self.initialize_config()

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.debug_level = config['debug_level']
        self.scale_histograms = config['scale_histograms']

        # Find pt-hat bin index
        self.pt_hat_bins = config['pt_hat_bins']
        self.n_pt_hat_bins = len(self.pt_hat_bins) - 1
        for i,bin in enumerate(self.pt_hat_bins):
            if bin == self.pt_hat_min:
                self.pt_hat_bin = i

        # Create output dir
        self.output_dir = os.path.join(self.output_dir, str(self.pt_hat_bin))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    # ---------------------------------------------------------------
    # Main processing function
    # ---------------------------------------------------------------
    def analyze_jetscape_events(self):

        print('Analyzing pt-hat: {} - {} ...'.format(self.pt_hat_min, self.pt_hat_max))

        # Read JETSCAPE output, get hadrons, do jet finding, and write histograms to ROOT file
        self.run_jetscape_analysis()

        # Scale histograms according to pthard bins cross-section
        # Note: if merging multiple files for given pt-hat bin, likely want to scale only at that point
        #       in order to scale by xsec/n_event
        if self.scale_histograms:
            print("Scaling pt-hat bin...")
            scale_histograms.scale_histograms(self.output_dir, self.pt_hat_bin, self.n_event_max, bRemoveOutliers=False)

    # ---------------------------------------------------------------
    # Main processing function for a single pt-hat bin
    # ---------------------------------------------------------------
    def run_jetscape_analysis(self):

        # Initialize output objects
        self.initialize_output_objects()

        # Read chunk of events into a dataframe
        # Fields: particle_ID, status, E, px, py, pz
        df_event_chunk = pd.read_parquet(self.input_file_hadrons)
        self.n_event_max = df_event_chunk.shape[0]

        # Iterate through events
        self.analyze_event_chunk(df_event_chunk)

        # Write analysis task output to ROOT file
        self.write_output_objects()

    # ---------------------------------------------------------------
    # Analyze event chunk
    # ---------------------------------------------------------------
    def analyze_event_chunk(self, df_event_chunk):

        # Loop through events
        for i,event in df_event_chunk.iterrows():

            if i % 1000 == 0:
                print('event: {}'.format(i))

            # Call user-defined function to analyze event
            self.analyze_event(event)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_output_objects(self):

        # Event histograms
        self.hNevents = ROOT.TH1F('hNevents', 'hNevents', self.n_pt_hat_bins, 0, self.n_pt_hat_bins)
        self.hCrossSection = ROOT.TH1F('hCrossSection', 'hCrossSection', self.n_pt_hat_bins, 0, self.n_pt_hat_bins)

        # Initialize user-defined output objects
        self.initialize_user_output_objects()

    # ---------------------------------------------------------------
    # Save all ROOT histograms and trees to file
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Fill cross-section
        self.hCrossSection.SetBinContent(self.pt_hat_bin+1, self.pt_hat_xsec)
        self.hCrossSection.SetBinError(self.pt_hat_bin+1, self.pt_hat_xsec_err)

        # Set N events
        self.hNevents.SetBinContent(self.pt_hat_bin+1, self.n_event_max)

        # Save output objects
        outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
        fout = ROOT.TFile(outputfilename, 'recreate')
        fout.cd()
        for attr in dir(self):

            obj = getattr(self, attr)

            # Write all ROOT histograms and trees to file
            types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
            if isinstance(obj, types):
                obj.Write()
                obj.SetDirectory(0)
                del obj

        fout.Close()

    # ---------------------------------------------------------------
    # Fill hadrons into vector of fastjet pseudojets
    #
    # By default, select all particles
    # If select_status='+', select only positive status particles
    # If select_status='-', select only positive status particles
    # ---------------------------------------------------------------
    def fill_fastjet_constituents(self, event, select_status=None, select_charged=False):
    
        if select_status == '-':
            status_mask = (event['status'] < 0)
        elif select_status == '+':
            status_mask = (event['status'] > -1)
        else:
            # Picked a value to make an all true mask. We don't select anything
            status_mask = event["status"] > -1e6

        # Default to an all true mask
        charged_mask = np.ones(len(status_mask)) > 0
        if select_charged:
            # NOTE: This is super inefficient - we can seriously accelerate this with numba.
            #       But this apparently works for now, so we leave it as is for now.
            # Create an all false mask. We'll fill it in with the charged constituents
            charged_mask = np.ones(len(event['particle_ID'])) < 0
            for i, pid_value in enumerate(event['particle_ID']):
                # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if np.abs(pid_value) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    charged_mask[i] = True

        full_mask = status_mask & charged_mask
        px = event['px'][full_mask]
        py = event['py'][full_mask]
        pz = event['pz'][full_mask]
        e = event['E'][full_mask]
        pid = event['particle_ID'][full_mask]

        # Create a vector of fastjet::PseudoJets from arrays of px,py,pz,e
        fj_particles = fjext.vectorize_px_py_pz_e(px, py, pz, e)

        # Set pid as user_index
        for i,p in enumerate(fj_particles):
            fj_particles[i].set_user_index(int(pid[i]))

        return fj_particles

    # ---------------------------------------------------------------
    # This function is called once per setting
    # You must implement this
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):
        raise NotImplementedError('You must implement initialize_user_output_objects()!')

    # ---------------------------------------------------------------
    # This function is called once per event (per setting)
    # You must implement this
    # ---------------------------------------------------------------
    def analyze_event(self, event):
        raise NotImplementedError('You must implement analyze_event()!')
