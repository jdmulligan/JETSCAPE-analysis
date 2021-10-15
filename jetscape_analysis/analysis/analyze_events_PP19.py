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
import yaml

# Analysis
import ROOT
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj

sys.path.append('../..')
from jetscape_analysis.analysis import analyze_events_base

################################################################
class AnalyzeJetscapeEvents_Example(analyze_events_base.AnalyzeJetscapeEvents_Base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_Example, self).__init__(config_file=config_file,
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

        # Hadron histograms
        hname = 'hChHadronPt_eta'
        pt_bins = [9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 140.0, 165.0, 250.0, 400.0]
        n_pt_bins = len(pt_bins) - 1
        pt_bin_array = array('d', pt_bins)
        eta_bins = [-5., -3., -1., 1., 3., 5.]
        n_eta_bins = len(eta_bins) - 1
        eta_bin_array = array('d', eta_bins)
        h = ROOT.TH2F(hname, hname, n_pt_bins, pt_bin_array, n_eta_bins, eta_bin_array)
        h.Sumw2()
        setattr(self, hname, h)
        
        hname = 'hD0Pt_eta'
        pt_bins = [2., 3., 4., 5., 6., 8., 10., 12.5, 15., 20., 25., 30., 40., 60., 100.]
        n_pt_bins = len(pt_bins) - 1
        pt_bin_array = array('d', pt_bins)
        eta_bins = [-5., -3., -1., 1., 3., 5.]
        n_eta_bins = len(eta_bins) - 1
        eta_bin_array = array('d', eta_bins)
        h = ROOT.TH2F(hname, hname, n_pt_bins, pt_bin_array, n_eta_bins, eta_bin_array)
        h.Sumw2()
        setattr(self, hname, h)
        
        hname = 'hHadronPID_eta{}'.format(self.abs_track_eta_max)
        setattr(self, hname, ROOT.TH1F(hname, hname, 10000, -5000, 5000))
        
        hname = 'hChHadronEtaPhi'
        setattr(self, hname, ROOT.TH2F(hname, hname, 100, -5, 5, 100, -3.2, 3.2))
        
        # Parton histograms
        hname = 'hPartonPt_eta{}'.format(self.abs_track_eta_max)
        setattr(self, hname, ROOT.TH1F(hname, hname, 300, 0.0, 300.0))
        
        hname = 'hPartonPID_eta{}'.format(self.abs_track_eta_max)
        setattr(self, hname, ROOT.TH1F(hname, hname, 10000, -5000, 5000))
        
        hname = 'hPartonEtaPhi'
        setattr(self, hname, ROOT.TH2F(hname, hname, 100, -5, 5, 100, -3.2, 3.2))

        # Jet histograms
        for jetR in self.jetR_list:

            hname = 'hJetPt_eta_R{}'.format(jetR)
            h = ROOT.TH2F(hname, hname, 1000, 0, 1000, 60, -3.0, 3.0)
            setattr(self, hname, h)

            hname = 'hJetEtaPhi_R{}'.format(jetR)
            h = ROOT.TH2F(hname, hname, 100, -5, 5, 100, 0, 6.28)
            setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Get list of hadrons from the event, and fill some histograms
        hadrons = event.hadrons(min_track_pt=self.min_track_pt)
        self.fill_hadron_histograms(hadrons)

        # Get list of final-state partons from the event, and fill some histograms
        partons = event.final_partons()
        self.fill_parton_histograms(partons)

        # Create list of fastjet::PseudoJets
        fj_hadrons = []
        fj_hadrons = self.fill_fastjet_constituents(hadrons)

        # Loop through specified jet R
        for jetR in self.jetR_list:

            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(5.)
            if self.debug_level > 0:
                print('jet definition is:', jet_def)
                print('jet selector is:', jet_selector, '\n')

            # Do jet finding
            jets = []
            jets_selected = []
            cs = fj.ClusterSequence(fj_hadrons, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # Fill some jet histograms
            self.fill_jet_histograms(jets_selected, jetR)

    # ---------------------------------------------------------------
    # Fill hadron histograms
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, hadrons):
    
        # Loop through hadrons
        for hadron in hadrons:

            # Fill some basic hadron info
            pid = hadron.pid

            momentum = hadron.momentum
            pt = momentum.pt()
            eta = momentum.eta()
            phi = momentum.phi()  # [-pi, pi]

            # Fill charged particle histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            # (assuming weak strange decays are off, but charm decays are on)
            # Neglect e+, mu+ (11, 13)
            if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                getattr(self, 'hChHadronPt_eta').Fill(pt, eta, 1/pt) # Fill with weight 1/pt, to form 1/pt dN/dpt
                getattr(self, 'hChHadronEtaPhi').Fill(eta, phi)

            if abs(eta) < self.abs_track_eta_max:
            
                hname = 'hHadronPID_eta{}'.format(self.abs_track_eta_max)
                getattr(self, hname).Fill(pid)
                
            if abs(pid) == 421:
            
                getattr(self, 'hD0Pt_eta').Fill(pt, eta, 1/pt) # Fill with weight 1/pt, to form 1/pt dN/dpt

    # ---------------------------------------------------------------
    # Fill final-state parton histograms
    # ---------------------------------------------------------------
    def fill_parton_histograms(self, partons):

        # Loop through partons
        for parton in partons:

            # Fill some basic parton info
            pid = parton.pid

            momentum = parton.momentum
            pt = momentum.pt()
            eta = momentum.eta()
            phi = momentum.phi()  # [-pi, pi]
            
            getattr(self, 'hPartonEtaPhi').Fill(eta, phi)

            if abs(eta) < self.abs_track_eta_max:
            
                hname = 'hPartonPt_eta{}'.format(self.abs_track_eta_max)
                getattr(self, hname).Fill(pt)
                
                hname = 'hPartonPID_eta{}'.format(self.abs_track_eta_max)
                getattr(self, hname).Fill(pid)
    
    # ---------------------------------------------------------------
    # Fill jet histograms
    # ---------------------------------------------------------------
    def fill_jet_histograms(self, jets, jetR):

        for jet in jets:

            jet_pt = jet.pt()
            jet_eta = jet.eta()
            jet_phi = jet.phi()  # [0, 2pi]

            getattr(self, 'hJetPt_eta_R{}'.format(jetR)).Fill(jet_pt, jet_eta)
            getattr(self, 'hJetEtaPhi_R{}'.format(jetR)).Fill(jet_eta, jet_phi)

##################################################################
if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser(description="Generate JETSCAPE events")
    parser.add_argument(
        "-c",
        "--configFile",
        action="store",
        type=str,
        metavar="configFile",
        default="/home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml",
        help="Path of config file for analysis",
    )
    parser.add_argument(
        "-i",
        "--inputDir",
        action="store",
        type=str,
        metavar="inputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Input directory containing JETSCAPE output files",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        type=str,
        metavar="outputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Output directory for output to be written to",
    )

    # Parse the arguments
    args = parser.parse_args()

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File "{0}" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid inputDir is given, exit
    if not os.path.exists(args.inputDir):
        print('File "{0}" does not exist! Exiting!'.format(args.inputDir))
        sys.exit(0)

    analysis = AnalyzeJetscapeEvents_Example(config_file=args.configFile, input_dir=args.inputDir, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
