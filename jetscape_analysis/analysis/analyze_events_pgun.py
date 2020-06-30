#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE output file

  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# General
import sys
import os
import argparse

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjext
import ROOT
import yaml
from array import *

sys.path.append('../..')
sys.path.append('.')
from jetscape_analysis.analysis import analyze_events
from jetscape_analysis.analysis.event import event_hepmc
from jetscape_analysis.base import common_base

################################################################
class AnalyzeJetscapeEvents_PGun(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_PGun, self).__init__(config_file=config_file,
                                                            input_file=input_file,
                                                            output_dir=output_dir,
                                                            **kwargs)
        self.initialize_user_config()
        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_user_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)
        
        self.min_track_pt = config['min_track_pt']
        self.abs_track_eta_max = config['abs_track_eta_max']
        
        self.jetR_list = config['jetR']
        self.min_jet_pt = config['min_jet_pt']

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):

        # Final-state parton histograms
        self.hPartonN = ROOT.TH1F("hPartonN", "hPartonN", 100, 0, 100)
        self.hPartonPt = ROOT.TH1F("hPartonPt", "hPartonPt", 100, 0.0, 100.0)
        self.hPartonE = ROOT.TH1F("hPartonE", "hPartonE", 100, 0.0, 100.0)
        self.hPartonPID = ROOT.TH1F("hPartonPID", "hPartonPID", 10000, -5000, 5000)
        self.hPartonEtaPhi = ROOT.TH2F("hPartonEtaPhi", "hPartonEtaPhi", 100, -5, 5, 100, -3.2, 3.2)

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Get list of final-state partons from the event, and fill some histograms
        partons = event.final_partons()
        self.fill_parton_histograms(partons)

    # ---------------------------------------------------------------
    # Fill final-state parton histograms
    # ---------------------------------------------------------------
    def fill_parton_histograms(self, partons):

        # Loop through partons
        for parton in partons:

            # Fill some basic parton info
            pid = parton.pid

            momentum = parton.momentum
            e = momentum.e
            pt = momentum.pt()
            eta = momentum.eta()
            phi = momentum.phi()  # [-pi, pi]

            if pid == 21:
                self.hPartonE.Fill(e)
                self.hPartonPt.Fill(pt)
                self.hPartonEtaPhi.Fill(eta, phi)
            self.hPartonPID.Fill(pid)

        self.hPartonN.Fill(len(partons))
