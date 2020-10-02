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
import numpy as np

# Fastjet via python (from external library heppy)
import fastjet as fj
import ROOT

sys.path.append('../..')
from jetscape_analysis.analysis import analyze_events_base_PHYS

################################################################
class AnalyzeJetscapeEvents_PHYS(analyze_events_base_PHYS.AnalyzeJetscapeEvents_BasePHYS):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_PHYS, self).__init__(config_file=config_file,
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
        
        # Charged particle parameters
        self.charged_particle_eta_cut = config['charged_particle_eta_cut']

        self.file_CMS_hadron_0_5 = config['CMS_hadron_0_5']
        self.file_CMS_hadron_5_10 = config['CMS_hadron_5_10']
        self.file_CMS_hadron_30_50 = config['CMS_hadron_30_50']
        self.file_ATLAS_hadron = config['ATLAS_hadron']
        self.file_ALICE_hadron = config['ALICE_hadron']
        
        # Get binnings from data
        f = ROOT.TFile(self.file_CMS_hadron_0_5, 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.bins_CMS_hadron = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        f = ROOT.TFile(self.file_ALICE_hadron, 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.bins_ALICE_hadron = np.array(h.GetXaxis().GetXbins())
        f.Close()

        # Jet parameters
        self.min_track_pt = config['min_track_pt']
        self.jetR_list = config['jetR']
        self.min_jet_pt = config['min_jet_pt']

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):

        # Hadron histograms
        hname = 'hChargedPt_CMS'
        h = ROOT.TH1F(hname, hname, len(self.bins_CMS_hadron)-1, self.bins_CMS_hadron)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hChargedPt_ATLAS'
        h = ROOT.TH1F(hname, hname, 3000, 0, 300)
        h.Sumw2()
        setattr(self, hname, h)
        
        hname = 'hChargedPt_ALICE'
        h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_hadron)-1, self.bins_ALICE_hadron)
        h.Sumw2()
        setattr(self, hname, h)

        # Jet histograms
        for jetR in self.jetR_list:

            hname = 'hJetPt_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, 300, 0, 300)
            setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):
    
        # Get list of hadrons from the event, and fill some histograms
        hadrons = event.hadrons(min_track_pt=self.min_track_pt)
        self.fill_hadron_histograms(hadrons)

        # Create list of fastjet::PseudoJets
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
            cs = fj.ClusterSequence(fj_hadrons, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # Fill some jet histograms
            self.fill_jet_histograms(jets_selected, jetR)

    # ---------------------------------------------------------------
    # Fill hadron histograms
    # (assuming weak strange decays are off, but charm decays are on)
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, hadrons):
    
        # Loop through hadrons
        for hadron in hadrons:

            # Fill some basic hadron info
            pid = hadron.pid
            pt = hadron.momentum.pt()
            eta = hadron.momentum.eta()

            # CMS
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.charged_particle_eta_cut[0]:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_CMS').Fill(pt)

            # ATLAS
            # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            # Exclude e+, mu+ (11, 13)
            if abs(eta) < self.charged_particle_eta_cut[1]:
                if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_ATLAS').Fill(pt)
                    
            # ALICE
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.charged_particle_eta_cut[2]:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_ALICE').Fill(pt)
    
    # ---------------------------------------------------------------
    # Fill jet histograms
    # ---------------------------------------------------------------
    def fill_jet_histograms(self, jets, jetR):

        for jet in jets:
            jet_pt = jet.pt()
            getattr(self, 'hJetPt_R{}'.format(jetR)).Fill(jet_pt)

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
        "--inputFile",
        action="store",
        type=str,
        metavar="inputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/test.out",
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
    if not os.path.exists(args.inputFile):
        print('File "{0}" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    analysis = AnalyzeJetscapeEvents_PHYS(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
