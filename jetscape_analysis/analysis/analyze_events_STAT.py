#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE parquet output file,
  and write out a new parquet file containing calculated observables

  Author: James Mulligan (james.mulligan@berkeley.edu)
  Author: Raymond Ehlers (raymond.ehlers@cern.ch)
  """

from __future__ import print_function

# General
import itertools
import sys
import os
import argparse
import yaml
import numpy as np
import pandas as pd

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

sys.path.append('.')
from jetscape_analysis.analysis import analyze_events_base_STAT

################################################################
class AnalyzeJetscapeEvents_STAT(analyze_events_base_STAT.AnalyzeJetscapeEvents_BaseSTAT):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_STAT, self).__init__(config_file=config_file,
                                                         input_file=input_file,
                                                         output_dir=output_dir,
                                                         **kwargs)
        # Initialize config file
        self.initialize_user_config()
        
        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_user_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)
            
        self.sqrts = config['sqrt_s']
        self.output_file = config['output_file']
         
        # Load observable blocks
        self.hadron_observables = config['hadron']
        self.hadron_correlation_observables = config['hadron_correlations']
        self.inclusive_chjet_observables = config['inclusive_chjet']
        self.inclusive_jet_observables = None
        self.semi_inclusive_chjet_observables = None
        self.dijet_observables = None
        if 'inclusive_jet' in config:
            self.inclusive_jet_observables = config['inclusive_jet']
        if 'semi_inclusive_chjet' in config:
            self.semi_inclusive_chjet_observables = config['semi_inclusive_chjet']
        if 'dijet' in config:
            self.dijet_observables = config['dijet']
        
        # General jet finding parameters
        self.jet_R = config['jet_R']
        self.min_jet_pt = config['min_jet_pt']
        self.max_jet_y = config['max_jet_y']
        
        # General grooming parameters
        if 'SoftDrop' in config:
            self.grooming_settings = config['SoftDrop']
                                
    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    #
    # The jet finding is done on positive status particles (shower+recoil),
    # and the negative status particles (holes) are then used after jet
    # finding to perform corrections
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Create list of fastjet::PseudoJets (separately for jet shower particles and holes)
        fj_hadrons_positive = self.fill_fastjet_constituents(event, select_status='+')
        fj_hadrons_negative = self.fill_fastjet_constituents(event, select_status='-')
        
        # Create list of charged particles
        fj_hadrons_positive_charged = self.fill_fastjet_constituents(event, select_status='+',
                                                                     select_charged=True)
        fj_hadrons_negative_charged = self.fill_fastjet_constituents(event, select_status='-',
                                                                     select_charged=True)
        
        # Fill hadron histograms for jet shower particles
        self.fill_hadron_histograms(fj_hadrons_positive, status='+')

        # Loop through specified jet R
        for jetR in self.jet_R:
        
            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(self.max_jet_y)

            # Full jets
            # -----------------
            cs = fj.ClusterSequence(fj_hadrons_positive, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # Fill inclusive full jet histograms
            [self.analyze_inclusive_jet(jet, fj_hadrons_positive, fj_hadrons_negative, jetR, charged=False) for jet in jets_selected]
            
            # Charged jets
            # -----------------
            cs_charged = fj.ClusterSequence(fj_hadrons_positive_charged, jet_def)
            jets_charged = fj.sorted_by_pt(cs_charged.inclusive_jets())
            jets_selected_charged = jet_selector(jets_charged)

            # Fill inclusive charged jet histograms
            [self.analyze_inclusive_jet(jet, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR, charged=True) for jet in jets_selected_charged]
            
            # Fill semi-inclusive jet correlations
            if self.semi_inclusive_chjet_observables:
                self.fill_semi_inclusive_chjet_histograms(jets_selected_charged, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR)
            
            # Fill dijet histograms
            if self.dijet_observables:
                self.fill_dijet_histograms(jets_selected_charged, fj_hadrons_negative_charged, jetR)
            
        # Fill the observables dict to a new entry in the event list
        self.output_event_list.append(self.observable_dict_event)

    # ---------------------------------------------------------------
    # Fill hadron histograms
    # (assuming weak strange decays are off, but charm decays are on)
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, fj_particles, status='+'):
    
        # Create empty lists of hadrons
        if self.sqrts in [2760, 5020]:
            hadron_pt_ch_alice = []
            hadron_pt_pi_alice = []
            hadron_pt_pi0_alice = []
            hadron_pt_ch_atlas = []
            hadron_pt_ch_cms = []
        elif self.sqrts in [200]:
            hadron_pt_pi0_phenix = []
            hadron_pt_ch_star = []

        # Loop through hadrons
        for i,particle in enumerate(fj_particles):

            # Fill some basic hadron info
            pid = particle.user_index()
            pt = particle.pt()
            eta = particle.eta()

            if self.sqrts in [2760, 5020]:

                # ALICE
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                pt_min = self.hadron_observables['pt_ch_alice']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_alice']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_alice']['eta_cut']:
                        if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                            hadron_pt_ch_alice.append(pt)
                        
                # Charged pion
                pt_min = self.hadron_observables['pt_pi_alice']['pt'][0]
                pt_max = self.hadron_observables['pt_pi_alice']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_pi_alice']['eta_cut']:
                        if abs(pid) == 211:
                            hadron_pt_pi_alice.append(pt)
                          
                # Neutral pions
                if self.sqrts in [2760]:
                    pt_min = self.hadron_observables['pt_pi0_alice']['pt'][0]
                    pt_max = self.hadron_observables['pt_pi0_alice']['pt'][1]
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_pi0_alice']['eta_cut']:
                            if abs(pid) == 111:
                                hadron_pt_pi0_alice.append(pt)

                # ATLAS
                # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                # Exclude e+, mu+ (11, 13)
                pt_min = self.hadron_observables['pt_ch_atlas']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_atlas']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_atlas']['eta_cut']:
                        if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                            hadron_pt_ch_atlas.append(pt)

                # CMS
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                pt_min = self.hadron_observables['pt_ch_cms']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_cms']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_cms']['eta_cut']:
                        if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                            hadron_pt_ch_cms.append(pt)
                            
            elif self.sqrts in [200]:
            
                # PHENIX
                # Neutral pions
                pt_min = self.hadron_observables['pt_pi0_phenix']['pt'][0]
                pt_max = self.hadron_observables['pt_pi0_phenix']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_pi0_phenix']['eta_cut']:
                        if abs(pid) == 111:
                            hadron_pt_pi0_phenix.append(pt)
            
                # STAR
                # Charged hadrons (pi+, K+, p+)
                pt_min = self.hadron_observables['pt_ch_star']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_star']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_star']['eta_cut']:
                        if abs(pid) in [211, 321, 2212]:
                            hadron_pt_ch_star.append(pt)

        # Fill lists to output dictionary
        if self.sqrts in [2760, 5020]:
            self.observable_dict_event['hadron_pt_ch_alice'] = hadron_pt_ch_alice
            self.observable_dict_event['hadron_pt_pi_alice'] = hadron_pt_pi_alice
            if self.sqrts in [2760]:
                self.observable_dict_event['hadron_pt_pi0_alice'] = hadron_pt_pi0_alice
            self.observable_dict_event['hadron_pt_ch_atlas'] = hadron_pt_ch_atlas
            self.observable_dict_event['hadron_pt_ch_cms'] = hadron_pt_ch_cms
        elif self.sqrts in [200]:
            self.observable_dict_event['hadron_pt_pi0_phenix'] = hadron_pt_pi0_phenix
            self.observable_dict_event['hadron_pt_ch_star'] = hadron_pt_ch_star

    # ---------------------------------------------------------------
    # Fill inclusive jet histograms
    #
    # To correct jet pt: sum up the hole pt within R of jet axis
    # To correct substructure:
    #   - For additive observables, sum up the hole substructure observable and subtract
    #   - For identified objects within jets or between jets (groomed jet, subjet, jet axis, delta_phi),
    #     construct the observable only from the positive status particles, and correct only the jet pt
    # ---------------------------------------------------------------
    def analyze_inclusive_jet(self, jet, fj_hadrons_positive, fj_hadrons_negative, jetR, charged=False):

        # Get the corrected jet pt by subtracting the negative recoils within R
        holes_in_jet = []
        negative_pt = 0.
        for hadron in fj_hadrons_negative:
            if jet.delta_R(hadron) < jetR:
                negative_pt += hadron.pt()
                holes_in_jet.append(hadron)
    
        jet_pt_uncorrected = jet.pt()               # uncorrected pt: shower+recoil
        jet_pt = jet_pt_uncorrected - negative_pt   # corrected pt: shower+recoil-holes
        
        # Construct groomed jet
        gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
        jet_groomed_lund = gshop.soft_drop(self.grooming_settings[0]['beta'], self.grooming_settings[0]['zcut'], jetR)

        # Fill histograms
        if charged:
            self.fill_charged_jet_histograms(jet, jet_groomed_lund, holes_in_jet, jet_pt, jetR)
        else:
            self.fill_full_jet_histograms(jet, jet_pt, jetR)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_full_jet_histograms(self, jet, jet_pt, jetR):

        # ALICE RAA
        if abs(jet.eta()) < (self.inclusive_jet_observables['pt_alice']['eta_cut_R'] - jetR):

            # Check leading track requirement
            if jetR == 0.2:
                min_leading_track_pt = 5.
            else:
                min_leading_track_pt = 7.

            accept_jet = False
            for constituent in jet.constituents():
                if constituent.pt() > min_leading_track_pt:
                    # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                    if abs(constituent.user_index()) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                        accept_jet = True

            if accept_jet:
                self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}'] = jet_pt
            self.observable_dict_event[f'inclusive_jet_pt_alice_no_ptlead_cut_R{jetR}'] = jet_pt

        # ATLAS RAA
        if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
            self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}'] = jet_pt

        # CMS RAA
        if abs(jet.eta()) < self.inclusive_jet_observables['pt_cms']['eta_cut']:
            self.observable_dict_event[f'inclusive_jet_pt_cms_R{jetR}'] = jet_pt
    
    # ---------------------------------------------------------------
    # Fill inclusive charged jet histograms
    # ---------------------------------------------------------------
    def fill_charged_jet_histograms(self, jet, jet_groomed_lund, holes_in_jet, jet_pt, jetR):
    
        if self.sqrts == 2760:

            # g
            if 40 < jet_pt < 60 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
                g = 0
                for constituent in jet.constituents():
                    g += constituent.pt() / jet_pt * constituent.delta_R(jet)
                for hadron in holes_in_jet:
                    g -= hadron.pt() / jet_pt * hadron.delta_R(jet)
                self.observable_dict_event[f'inclusive_chjet_g_alice_R{jetR}'] = g

            # Jet mass
            if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            
                jet_mass = jet.m()

                # Add holes together as four vectors, and then calculate their mass, removing it from the jet mass.
                if holes_in_jet:
                    # Avoid modifying the original hole.
                    hole_four_vector = fj.PseudoJet()
                    for hadron in holes_in_jet:
                        hole_four_vector += hadron

                    # Remove mass from holes
                    jet_mass -= hole_four_vector.m()

                self.observable_dict_event[f'inclusive_chjet_mass_alice_R{jetR}'] = jet_mass

        elif self.sqrts == 5020:

            # Soft Drop
            if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['zg_alice']['eta_cut_R'] - jetR):
                if jet_groomed_lund:
                    theta_g = jet_groomed_lund.Delta() / jetR
                    zg = jet_groomed_lund.z()
                    # Note: untagged jets will return negative value
                self.observable_dict_event[f'inclusive_chjet_zg_alice_R{jetR}'] = zg
                self.observable_dict_event[f'inclusive_chjet_tg_alice_R{jetR}'] = theta_g

    # ---------------------------------------------------------------
    # Fill semi-inclusive charged jet histograms
    #
    # Note: We may need a lower jet pt range to determine cref, but I didn't look into this.
    # Note: Doesn't account for detector effects on hadron.
    # ---------------------------------------------------------------
    def fill_semi_inclusive_chjet_histograms(self, jets_selected, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR):

        # Define trigger classes for both traditional h-jet analysis and Nsubjettiness analysis
        hjet_low_trigger_range_276 = self.semi_inclusive_chjet_observables['hjet_alice']['low_trigger_range_276']
        hjet_low_trigger_range_502 = self.semi_inclusive_chjet_observables['hjet_alice']['low_trigger_range_502']
        hjet_high_trigger_range = self.semi_inclusive_chjet_observables['hjet_alice']['high_trigger_range']
        nsubjettiness_low_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['low_trigger_range']
        nsubjettiness_high_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['high_trigger_range']
        
        # Define Nsubjettiness calculators
        axis_definition = fjcontrib.KT_Axes()
        measure_definition = fjcontrib.UnnormalizedMeasure(1)
        n_subjettiness_calculator1 = fjcontrib.Nsubjettiness(1, axis_definition, measure_definition)
        n_subjettiness_calculator2 = fjcontrib.Nsubjettiness(2, axis_definition, measure_definition)

        for hadron in fj_hadrons_positive_charged:
        
            if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['hjet_alice']['hadron_eta_cut']:

                # Search for hadron trigger
                hjet_found_low_276 = False
                hjet_found_low_502 = False
                hjet_found_high = False
                nsubjettiness_found_low = False
                nsubjettiness_found_high = False

                if hjet_low_trigger_range_276[0] < hadron.pt() < hjet_low_trigger_range_276[1]:
                    hjet_found_low_276 = True
                if hjet_low_trigger_range_502[0] < hadron.pt() < hjet_low_trigger_range_502[1]:
                    hjet_found_low_502 = True
                if hjet_high_trigger_range[0] < hadron.pt() < hjet_high_trigger_range[1]:
                    hjet_found_high = True
                if nsubjettiness_low_trigger_range[0] < hadron.pt() < nsubjettiness_low_trigger_range[1]:
                    nsubjettiness_found_low = True
                if nsubjettiness_high_trigger_range[0] < hadron.pt() < nsubjettiness_high_trigger_range[1]:
                    nsubjettiness_found_high = True
                found_trigger =  hjet_found_low_276 or hjet_found_low_502 or hjet_found_high or nsubjettiness_found_low or nsubjettiness_found_high

                # Record N triggers
                self.observable_dict_event[f'semi_inclusive_chjet_hjet_ntrigger_alice_R{jetR}'] = hadron.pt()
                self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_ntrigger_alice_R{jetR}'] = hadron.pt()
                
                # Search for recoil jets
                if found_trigger:
                    for jet in jets_selected:
                        if abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
                                        
                            # Get the corrected jet pt: shower+recoil-holes
                            jet_pt = jet.pt()
                            for temp_hadron in fj_hadrons_negative_charged:
                                if jet.delta_R(temp_hadron) < jetR:
                                    jet_pt -= temp_hadron.pt()

                            # Jet yield and Delta phi
                            if hjet_found_low_276:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    self.observable_dict_event[f'semi_inclusive_chjet_IAA_lowTrigger_alice_R{jetR}_276'] = jet_pt

                                if 40 < jet_pt < 60:
                                    self.observable_dict_event[f'semi_inclusive_chjet_dphi_lowTrigger_alice_R{jetR}_276'] = np.abs(hadron.delta_phi_to(jet))
                                    
                            if hjet_found_high:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    self.observable_dict_event[f'semi_inclusive_chjet_IAA_highTrigger_alice_R{jetR}_276'] = jet_pt

                                if 40 < jet_pt < 60:
                                    self.observable_dict_event[f'semi_inclusive_chjet_dphi_highTrigger_alice_R{jetR}_276'] = np.abs(hadron.delta_phi_to(jet))

                            # Nsubjettiness
                            if nsubjettiness_found_low:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    if 40 < jet_pt < 60:
                                        tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                        tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                        if tau1 > 1e-3:
                                            self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_lowTrigger_alice_R{jetR}'] = tau2/tau1
                                            
                            if nsubjettiness_found_high:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    if 40 < jet_pt < 60:
                                        tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                        tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                        if tau1 > 1e-3:
                                            self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_highTrigger_alice_R{jetR}'] = tau2/tau1

    #---------------------------------------------------------------
    # Return leading jet (or subjet)
    #---------------------------------------------------------------
    def leading_jet(self, jets):

        leading_jet = None
        for jet in jets:

            if not leading_jet:
                leading_jet = jet
            
            if jet.pt() > leading_jet.pt():
                leading_jet = jet

        return leading_jet

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

    analysis = AnalyzeJetscapeEvents_STAT(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
