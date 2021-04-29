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
from pathlib import Path

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
        # Update the output_file to contain the labeling in the final_state_hadrons file.
        # We use this naming convention as the flag for whether we should attempt to rename it.
        if "final_state_hadrons" in self.input_file_hadrons:
            _input_filename = Path(self.input_file_hadrons).name
            # The filename will be something like "observables_{sqrts}_0000_00.parquet", assuming
            # that the original name was "observables_{sqrts}"
            self.output_file = _input_filename.replace("final_state_hadrons", self.output_file)
            print(f'Updated output_file name to "{self.output_file}" in order to add identifying indices.')

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
        else:
            self.grooming_settings = None
                                
    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    #
    # The jet finding is done on positive status particles (shower+recoil),
    # and the negative status particles (holes) are then used after jet
    # finding to perform corrections
    # ---------------------------------------------------------------
    def analyze_event(self, event):
    
        # Initialize empty list for each output observable
        self.initialize_output_lists()

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
        self.fill_hadron_histograms(fj_hadrons_negative, status='-')        

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
            [self.analyze_inclusive_jet(jet, fj_hadrons_positive, fj_hadrons_negative, jetR, full_jet=True) for jet in jets_selected]
            
            # Charged jets
            # -----------------
            cs_charged = fj.ClusterSequence(fj_hadrons_positive_charged, jet_def)
            jets_charged = fj.sorted_by_pt(cs_charged.inclusive_jets())
            jets_selected_charged = jet_selector(jets_charged)

            # Fill inclusive charged jet histograms
            [self.analyze_inclusive_jet(jet, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR, full_jet=False) for jet in jets_selected_charged]
            
            # Fill semi-inclusive jet correlations
            if self.semi_inclusive_chjet_observables:
                if self.sqrts == 2760:
                    jetR_list = self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']+self.semi_inclusive_chjet_observables['nsubjettiness_alice']['jet_R']
                elif self.sqrts == 200:
                    jetR_list = self.semi_inclusive_chjet_observables['IAA_star']['jet_R']
                if jetR in jetR_list:
                    self.fill_semi_inclusive_chjet_histograms(jets_selected_charged, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR)
            
            # Fill dijet histograms
            if self.dijet_observables:
                self.fill_dijet_histograms(jets_selected_charged, fj_hadrons_negative_charged, jetR)
            
        # Fill the observables dict to a new entry in the event list
        self.output_event_list.append(self.observable_dict_event)

    # ---------------------------------------------------------------
    # Initialize empty list for each output observable
    # ---------------------------------------------------------------
    def initialize_output_lists(self):
    
        for observable in self.hadron_observables:
            self.observable_dict_event[f'hadron_{observable}'] = []
            self.observable_dict_event[f'hadron_{observable}_holes'] = []
            
        for observable in self.hadron_correlation_observables:
            self.observable_dict_event[f'hadron_correlations_{observable}'] = []
        
        if self.inclusive_jet_observables:
            for key,dict in self.inclusive_jet_observables.items():
                for jetR in dict['jet_R']:
                    if 'SoftDrop' in dict:
                        for grooming_setting in dict['SoftDrop']:
                            zcut = grooming_setting['zcut']
                            beta = grooming_setting['beta']
                            self.observable_dict_event[f'inclusive_jet_{key}_R{jetR}_zcut{zcut}_beta{beta}'] = []
                    else:
                        if 'charge_cms' in key:
                            for kappa in dict['kappa']:
                                self.observable_dict_event[f'inclusive_jet_{key}_R{jetR}_k{kappa}'] = []
                        else:
                            self.observable_dict_event[f'inclusive_jet_{key}_R{jetR}'] = []
        
        for key,dict in self.inclusive_chjet_observables.items():
            for jetR in dict['jet_R']:
                if 'SoftDrop' in dict:
                    for grooming_setting in dict['SoftDrop']:
                        zcut = grooming_setting['zcut']
                        beta = grooming_setting['beta']
                        if key == 'tg_alice' and jetR == 0.2 and zcut == 0.4:
                            continue
                        self.observable_dict_event[f'inclusive_chjet_{key}_R{jetR}_zcut{zcut}_beta{beta}'] = []
                else:
                    self.observable_dict_event[f'inclusive_chjet_{key}_R{jetR}'] = []
                    
        if self.semi_inclusive_chjet_observables:
            for key,dict in self.semi_inclusive_chjet_observables.items():
                for jetR in dict['jet_R']:
                    if self.sqrts == 2760:
                        self.observable_dict_event[f'semi_inclusive_chjet_{key}_R{jetR}_lowTrigger'] = []
                        self.observable_dict_event[f'semi_inclusive_chjet_{key}_R{jetR}_highTrigger'] = []
                        self.observable_dict_event[f'semi_inclusive_chjet_alice_trigger_pt'] = []
                    elif self.sqrts == 200:
                        self.observable_dict_event[f'semi_inclusive_chjet_{key}_R{jetR}'] = []
                        self.observable_dict_event[f'semi_inclusive_chjet_star_trigger_pt'] = []
            
        if self.dijet_observables:
            for key,dict in self.dijet_observables.items():
                for jetR in dict['jet_R']:
                    self.observable_dict_event[f'dijet_{key}_R{jetR}'] = []
                                    
    # ---------------------------------------------------------------
    # Fill hadron histograms
    # (assuming weak strange decays are off, but charm decays are on)
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, fj_particles, status='+'):
    
        # Note that for identified particles, we store holes of the identified species
        suffix = ''
        if status == '-':
            suffix = '_holes'
    
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
                            self.observable_dict_event[f'hadron_pt_ch_alice{suffix}'].append(pt)
                        
                # Charged pion
                pt_min = self.hadron_observables['pt_pi_alice']['pt'][0]
                pt_max = self.hadron_observables['pt_pi_alice']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_pi_alice']['eta_cut']:
                        if abs(pid) == 211:
                            self.observable_dict_event[f'hadron_pt_pi_alice{suffix}'].append(pt)
                          
                # Neutral pions
                if self.sqrts in [2760]:
                    pt_min = self.hadron_observables['pt_pi0_alice']['pt'][0]
                    pt_max = self.hadron_observables['pt_pi0_alice']['pt'][1]
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_pi0_alice']['eta_cut']:
                            if abs(pid) == 111:
                                self.observable_dict_event[f'hadron_pt_pi0_alice{suffix}'].append(pt)

                # ATLAS
                # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                # Exclude e+, mu+ (11, 13)
                pt_min = self.hadron_observables['pt_ch_atlas']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_atlas']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_atlas']['eta_cut']:
                        if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                            self.observable_dict_event[f'hadron_pt_ch_atlas{suffix}'].append(pt)

                # CMS
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                pt_min = self.hadron_observables['pt_ch_cms']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_cms']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_cms']['eta_cut']:
                        if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                            self.observable_dict_event[f'hadron_pt_ch_cms{suffix}'].append(pt)
                            
            elif self.sqrts in [200]:
            
                # PHENIX
                # Neutral pions
                pt_min = self.hadron_observables['pt_pi0_phenix']['pt'][0]
                pt_max = self.hadron_observables['pt_pi0_phenix']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_pi0_phenix']['eta_cut']:
                        if abs(pid) == 111:
                            self.observable_dict_event[f'hadron_pt_pi0_phenix{suffix}'].append(pt)
            
                # STAR
                # Charged hadrons (pi+, K+, p+)
                pt_min = self.hadron_observables['pt_ch_star']['pt'][0]
                pt_max = self.hadron_observables['pt_ch_star']['pt'][1]
                if pt > pt_min and pt < pt_max:
                    if abs(eta) < self.hadron_observables['pt_ch_star']['eta_cut']:
                        if abs(pid) in [211, 321, 2212]:
                            self.observable_dict_event[f'hadron_pt_ch_star{suffix}'].append(pt)

    # ---------------------------------------------------------------
    # Fill inclusive jet histograms
    #
    # To correct jet pt: sum up the hole pt within R of jet axis
    # To correct substructure:
    #   - For additive observables, sum up the hole substructure observable and subtract
    #   - For identified objects within jets or between jets (groomed jet, subjet, jet axis, delta_phi),
    #     construct the observable only from the positive status particles, and correct only the jet pt
    # ---------------------------------------------------------------
    def analyze_inclusive_jet(self, jet, fj_hadrons_positive, fj_hadrons_negative, jetR, full_jet=True):

        # Get the corrected jet pt by subtracting the negative recoils within R
        holes_in_jet = []
        negative_pt = 0.
        for hadron in fj_hadrons_negative:
            if jet.delta_R(hadron) < jetR:
                negative_pt += hadron.pt()
                holes_in_jet.append(hadron)
    
        jet_pt_uncorrected = jet.pt()               # uncorrected pt: shower+recoil
        jet_pt = jet_pt_uncorrected - negative_pt   # corrected pt: shower+recoil-holes
        
        # Fill histograms
        if full_jet:
        
            # Ungroomed
            self.fill_full_jet_ungroomed_histograms(jet, fj_hadrons_positive, holes_in_jet, jet_pt, jetR)
            
            # Groomed
            if self.grooming_settings:
                for grooming_setting in self.grooming_settings:
                    self.fill_full_jet_groomed_histograms(grooming_setting, jet, holes_in_jet, jet_pt, jetR)

        else:
                  
            # Ungroomed
            self.fill_charged_jet_ungroomed_histograms(jet, holes_in_jet, jet_pt, jetR)
        
            # Groomed
            if self.grooming_settings:
                for grooming_setting in self.grooming_settings:
                    self.fill_charged_jet_groomed_histograms(grooming_setting, jet, holes_in_jet, jet_pt, jetR)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_full_jet_ungroomed_histograms(self, jet, fj_hadrons_positive, holes_in_jet, jet_pt, jetR):
    
        if self.sqrts in [2760, 5020]:

            # ALICE RAA
            pt_min = self.inclusive_jet_observables['pt_alice']['pt'][0]
            pt_max = self.inclusive_jet_observables['pt_alice']['pt'][1]
            if jetR in self.inclusive_jet_observables['pt_alice']['jet_R']:
                if abs(jet.eta()) < (self.inclusive_jet_observables['pt_alice']['eta_cut_R'] - jetR):
                    if jet_pt > pt_min and jet_pt < pt_max:

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
                            self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}'].append(jet_pt)

            # ATLAS RAA
            pt_min = self.inclusive_jet_observables['pt_atlas']['pt'][0]
            pt_max = self.inclusive_jet_observables['pt_atlas']['pt'][1]
            if jetR in self.inclusive_jet_observables['pt_atlas']['jet_R']:
                if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
                    if jet_pt > pt_min and jet_pt < pt_max:
                        self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}'].append(jet_pt)
            
            if self.sqrts in [2760]:
                pt_min = self.inclusive_jet_observables['pt_rcp_atlas']['pt'][0]
                pt_max = self.inclusive_jet_observables['pt_rcp_atlas']['pt'][1]
                if jetR in self.inclusive_jet_observables['pt_rcp_atlas']['jet_R']:
                    if abs(jet.eta()) < self.inclusive_jet_observables['pt_rcp_atlas']['eta_cut']:
                        if jet_pt > pt_min and jet_pt < pt_max:
                            self.observable_dict_event[f'inclusive_jet_pt_rcp_atlas_R{jetR}'] .append(jet_pt)

            # CMS RAA
            pt_min = self.inclusive_jet_observables['pt_cms']['pt'][0]
            pt_max = self.inclusive_jet_observables['pt_cms']['pt'][1]
            if jetR in self.inclusive_jet_observables['pt_cms']['jet_R']:
                if abs(jet.eta()) < self.inclusive_jet_observables['pt_cms']['eta_cut']:
                    if jet_pt > pt_min and jet_pt < pt_max:
                        self.observable_dict_event[f'inclusive_jet_pt_cms_R{jetR}'].append(jet_pt)
                        
            # ATLAS D(z)
            pt_min = self.inclusive_jet_observables['Dz_atlas']['pt'][0]
            pt_max = self.inclusive_jet_observables['Dz_atlas']['pt'][1]
            if jetR in self.inclusive_jet_observables['Dz_atlas']['jet_R']:
                if abs(jet.rap()) < self.inclusive_jet_observables['Dz_atlas']['y_cut']:
                    if jet_pt > pt_min and jet_pt < pt_max:
                        for hadron in fj_hadrons_positive:
                            # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                            pid = hadron.user_index()
                            if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                if jet.delta_R(hadron) < jetR:
                                    z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                    self.observable_dict_event[f'inclusive_jet_Dz_atlas_R{jetR}'].append([jet_pt, z])
                                    self.observable_dict_event[f'inclusive_jet_Dpt_atlas_R{jetR}'].append([jet_pt, hadron.pt()])
            
            # CMS D(z)
            if self.sqrts == 2760:
                pt_min = self.inclusive_jet_observables['Dz_cms']['pt'][0]
                pt_max = self.inclusive_jet_observables['Dz_cms']['pt'][1]
                eta_range = self.inclusive_jet_observables['Dz_cms']['eta_cut']
                if jetR in self.inclusive_jet_observables['Dz_cms']['jet_R']:
                    if eta_range[0] < abs(jet.eta()) < eta_range[1]:
                        if jet_pt > pt_min and jet_pt < pt_max:
                            for hadron in fj_hadrons_positive:
                                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                pid = hadron.user_index()
                                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                    if jet.delta_R(hadron) < jetR:
                                        z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                        self.observable_dict_event[f'inclusive_jet_Dz_cms_R{jetR}'].append([jet_pt, z])
                                        self.observable_dict_event[f'inclusive_jet_Dpt_cms_R{jetR}'].append([jet_pt, hadron.pt()])
                             
            # CMS jet charge
            if self.sqrts == 5020:
                pt_min = self.inclusive_jet_observables['charge_cms']['pt_min']
                if jetR in self.inclusive_jet_observables['charge_cms']['jet_R']:
                    if abs(jet.eta()) < self.inclusive_jet_observables['charge_cms']['eta_cut']:
                        if jet_pt > pt_min:
                            for kappa in self.inclusive_jet_observables['charge_cms']['kappa']:
                                sum = 0
                                for hadron in fj_hadrons_positive:
                                    if hadron.pt() > self.inclusive_jet_observables['charge_cms']['track_pt_min']:
                                        # Charged particles (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                        pid = hadron.user_index()
                                        if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                            if jet.delta_R(hadron) < jetR:
                                                sum += self.charge(pid) * np.power(hadron.pt(), kappa)
                                charge = sum / np.power(jet_pt, kappa)
                                self.observable_dict_event[f'inclusive_jet_charge_cms_R{jetR}_k{kappa}'].append(charge)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_full_jet_groomed_histograms(self, grooming_setting, jet, holes_in_jet, jet_pt, jetR):
    
        # Construct groomed jet
        gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
        zcut = grooming_setting['zcut']
        beta = grooming_setting['beta']
        jet_groomed_lund = gshop.soft_drop(beta, zcut, jetR)
        if not jet_groomed_lund:
            return
        
        if self.sqrts == 5020:

            # CMS m_g
            if grooming_setting in self.inclusive_jet_observables['mg_cms']['SoftDrop']:
                pt_min = self.inclusive_jet_observables['mg_cms']['pt'][0]
                pt_max = self.inclusive_jet_observables['mg_cms']['pt'][1]
                if jetR in self.inclusive_jet_observables['mg_cms']['jet_R']:
                    if abs(jet.eta()) < (self.inclusive_jet_observables['mg_cms']['eta_cut']):
                        if jet_pt > pt_min and jet_pt < pt_max:
                            if jet_groomed_lund.Delta() > self.inclusive_jet_observables['mg_cms']['dR']:
                                mg = jet_groomed_lund.pair().m()
                                # Note: untagged jets will return negative value
                                self.observable_dict_event[f'inclusive_jet_mg_cms_R{jetR}_zcut{zcut}_beta{beta}'].append([jet_pt, mg])
            
            # CMS z_g
            if grooming_setting in self.inclusive_jet_observables['zg_cms']['SoftDrop']:
                pt_min = self.inclusive_jet_observables['zg_cms']['pt'][0]
                pt_max = self.inclusive_jet_observables['zg_cms']['pt'][1]
                if jetR in self.inclusive_jet_observables['zg_cms']['jet_R']:
                    if abs(jet.eta()) < (self.inclusive_jet_observables['zg_cms']['eta_cut']):
                        if jet_pt > pt_min and jet_pt < pt_max:
                            if jet_groomed_lund.Delta() > self.inclusive_jet_observables['zg_cms']['dR']:
                                zg = jet_groomed_lund.z()
                                # Note: untagged jets will return negative value
                                self.observable_dict_event[f'inclusive_jet_zg_cms_R{jetR}_zcut{zcut}_beta{beta}'].append([jet_pt, zg])
            
    # ---------------------------------------------------------------
    # Fill inclusive charged jet histograms
    # ---------------------------------------------------------------
    def fill_charged_jet_ungroomed_histograms(self, jet, holes_in_jet, jet_pt, jetR):
    
        if self.sqrts == 2760:
        
            # ALICE charged jet RAA
            pt_min = self.inclusive_chjet_observables['pt_alice']['pt_ch'][0]
            pt_max = self.inclusive_chjet_observables['pt_alice']['pt_ch'][1]
            if jetR in self.inclusive_chjet_observables['pt_alice']['jet_R']:
                if abs(jet.eta()) < (self.inclusive_chjet_observables['pt_alice']['eta_cut']):
                    if jet_pt > pt_min and jet_pt < pt_max:

                        # Check leading track requirement
                        accept_jet = False
                        for constituent in jet.constituents():
                            if constituent.pt() > self.inclusive_chjet_observables['pt_alice']['leading_track_min_pt']:
                                # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                if abs(constituent.user_index()) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                    accept_jet = True
                        if accept_jet:
                            self.observable_dict_event[f'inclusive_chjet_pt_alice_R{jetR}'].append(jet_pt)

            # g
            pt_min = self.inclusive_chjet_observables['g_alice']['pt_ch'][0]
            pt_max = self.inclusive_chjet_observables['g_alice']['pt_ch'][1]
            if abs(jet.eta()) < (self.inclusive_chjet_observables['g_alice']['eta_cut_R'] - jetR):
                if jetR in self.inclusive_chjet_observables['g_alice']['jet_R']:
                    if jet_pt > pt_min and jet_pt < pt_max:
                        g = 0
                        for constituent in jet.constituents():
                            g += constituent.pt() / jet_pt * constituent.delta_R(jet)
                        for hadron in holes_in_jet:
                            g -= hadron.pt() / jet_pt * hadron.delta_R(jet)
                        self.observable_dict_event[f'inclusive_chjet_g_alice_R{jetR}'].append(g)
                    
            # pTD
            pt_min = self.inclusive_chjet_observables['ptd_alice']['pt_ch'][0]
            pt_max = self.inclusive_chjet_observables['ptd_alice']['pt_ch'][1]
            if abs(jet.eta()) < (self.inclusive_chjet_observables['ptd_alice']['eta_cut_R'] - jetR):
                if jetR in self.inclusive_chjet_observables['ptd_alice']['jet_R']:
                    if jet_pt > pt_min and jet_pt < pt_max:
                        sum = 0
                        for constituent in jet.constituents():
                            sum += np.power(constituent.pt(), 2)
                        for hadron in holes_in_jet:
                            sum -= np.power(hadron.pt(), 2)
                        self.observable_dict_event[f'inclusive_chjet_ptd_alice_R{jetR}'].append(np.sqrt(sum) / jet_pt)

            # Jet mass
            pt_min = self.inclusive_chjet_observables['mass_alice']['pt_ch'][0]
            pt_max = self.inclusive_chjet_observables['mass_alice']['pt_ch'][1]
            if abs(jet.eta()) < (self.inclusive_chjet_observables['mass_alice']['eta_cut_R'] - jetR):
                if jetR in self.inclusive_chjet_observables['mass_alice']['jet_R']:
                    if jet_pt > pt_min and jet_pt < pt_max:
            
                        jet_mass = jet.m()

                        # Add holes together as four vectors, and then calculate their mass, removing it from the jet mass.
                        if holes_in_jet:
                            # Avoid modifying the original hole.
                            hole_four_vector = fj.PseudoJet()
                            for hadron in holes_in_jet:
                                hole_four_vector += hadron

                            # Remove mass from holes
                            jet_mass -= hole_four_vector.m()

                        self.observable_dict_event[f'inclusive_chjet_mass_alice_R{jetR}'].append(jet_mass)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_charged_jet_groomed_histograms(self, grooming_setting, jet, holes_in_jet, jet_pt, jetR):
    
        # Construct groomed jet
        gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
        zcut = grooming_setting['zcut']
        beta = grooming_setting['beta']
        jet_groomed_lund = gshop.soft_drop(beta, zcut, jetR)
        if not jet_groomed_lund:
            return
        
        if self.sqrts == 5020:

            # Soft Drop
            if grooming_setting in self.inclusive_chjet_observables['zg_alice']['SoftDrop']:
                pt_min = self.inclusive_chjet_observables['zg_alice']['pt_ch'][0]
                pt_max = self.inclusive_chjet_observables['zg_alice']['pt_ch'][1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['zg_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['zg_alice']['jet_R']:
                        if jet_pt > pt_min and jet_pt < pt_max:
                            theta_g = jet_groomed_lund.Delta() / jetR
                            zg = jet_groomed_lund.z()
                            # Note: untagged jets will return negative value
                            self.observable_dict_event[f'inclusive_chjet_zg_alice_R{jetR}_zcut{zcut}_beta{beta}'].append(zg)
                            self.observable_dict_event[f'inclusive_chjet_tg_alice_R{jetR}_zcut{zcut}_beta{beta}'].append(theta_g)
            
    # ---------------------------------------------------------------
    # Fill semi-inclusive charged jet histograms
    # ---------------------------------------------------------------
    def fill_semi_inclusive_chjet_histograms(self, jets_selected, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR):

        if self.sqrts == 2760:
        
            # Define trigger classes for both traditional h-jet analysis and Nsubjettiness analysis
            hjet_low_trigger_range = self.semi_inclusive_chjet_observables['IAA_alice']['low_trigger_range']
            hjet_high_trigger_range = self.semi_inclusive_chjet_observables['IAA_alice']['high_trigger_range']
            nsubjettiness_low_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['low_trigger_range']
            nsubjettiness_high_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['high_trigger_range']
            
            pt_IAA = self.semi_inclusive_chjet_observables['IAA_alice']['pt_ch']
            pt_dphi = self.semi_inclusive_chjet_observables['dphi_alice']['pt_ch']
            pt_nsubjettiness = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['pt_ch']
            
            # Define Nsubjettiness calculators
            axis_definition = fjcontrib.KT_Axes()
            measure_definition = fjcontrib.UnnormalizedMeasure(1)
            n_subjettiness_calculator1 = fjcontrib.Nsubjettiness(1, axis_definition, measure_definition)
            n_subjettiness_calculator2 = fjcontrib.Nsubjettiness(2, axis_definition, measure_definition)

            for hadron in fj_hadrons_positive_charged:
            
                if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['IAA_alice']['hadron_eta_cut']:

                    # Search for hadron trigger
                    hjet_found_low = False
                    hjet_found_high = False
                    nsubjettiness_found_low = False
                    nsubjettiness_found_high = False

                    if hjet_low_trigger_range[0] < hadron.pt() < hjet_low_trigger_range[1]:
                        hjet_found_low = True
                    if hjet_high_trigger_range[0] < hadron.pt() < hjet_high_trigger_range[1]:
                        hjet_found_high = True
                    if nsubjettiness_low_trigger_range[0] < hadron.pt() < nsubjettiness_low_trigger_range[1]:
                        nsubjettiness_found_low = True
                    if nsubjettiness_high_trigger_range[0] < hadron.pt() < nsubjettiness_high_trigger_range[1]:
                        nsubjettiness_found_high = True
                    found_trigger =  hjet_found_low or hjet_found_high or nsubjettiness_found_low or nsubjettiness_found_high
            
                    if found_trigger:
                    
                        # Record hadron pt for trigger normalization
                        if jetR == min(self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']):
                            self.observable_dict_event[f'semi_inclusive_chjet_alice_trigger_pt'].append(hadron.pt())
                    
                        # Search for recoil jets
                        for jet in jets_selected:
                            if abs(jet.eta()) < (self.semi_inclusive_chjet_observables['IAA_alice']['eta_cut_R'] - jetR):
                                            
                                # Get the corrected jet pt: shower+recoil-holes
                                jet_pt = jet.pt()
                                for temp_hadron in fj_hadrons_negative_charged:
                                    if jet.delta_R(temp_hadron) < jetR:
                                        jet_pt -= temp_hadron.pt()

                                # Jet yield and Delta phi
                                if jetR in self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']:
                                    if hjet_found_low:
                                        if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                            if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_lowTrigger'].append(jet_pt)

                                        if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                            self.observable_dict_event[f'semi_inclusive_chjet_dphi_alice_R{jetR}_lowTrigger'].append(np.abs(hadron.delta_phi_to(jet)))
                                            
                                    if hjet_found_high:
                                        if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                            if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_highTrigger'].append(jet_pt)

                                        if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                            self.observable_dict_event[f'semi_inclusive_chjet_dphi_alice_R{jetR}_highTrigger'].append(np.abs(hadron.delta_phi_to(jet)))

                                # Nsubjettiness
                                if jetR in self.semi_inclusive_chjet_observables['nsubjettiness_alice']['jet_R']:
                                    if nsubjettiness_found_low:
                                        if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                            if pt_nsubjettiness[0] < jet_pt < pt_nsubjettiness[1]:
                                                tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                                tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                                if tau1 > 1e-3:
                                                    self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_alice_R{jetR}_lowTrigger'].append(tau2/tau1)
                                                    
                                    if nsubjettiness_found_high:
                                        if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                            if pt_nsubjettiness[0] < jet_pt < pt_nsubjettiness[1]:
                                                tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                                tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                                if tau1 > 1e-3:
                                                    self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_alice_R{jetR}_highTrigger'].append(tau2/tau1)

        if self.sqrts == 200:
        
            hjet_trigger_range = self.semi_inclusive_chjet_observables['IAA_star']['trigger_range']
            pt_IAA = self.semi_inclusive_chjet_observables['IAA_star']['pt_ch']
            pt_dphi = self.semi_inclusive_chjet_observables['dphi_star']['pt_ch']
            
            for hadron in fj_hadrons_positive_charged:
            
                if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['IAA_star']['hadron_eta_cut']:

                    # Search for hadron trigger
                    found_trigger = False
                    if hjet_trigger_range[0] < hadron.pt() < hjet_trigger_range[1]:
                        found_trigger = True
            
                    if found_trigger:
                    
                        # Record hadron pt for trigger normalization
                        if jetR == min(self.semi_inclusive_chjet_observables['IAA_star']['jet_R']):
                            self.observable_dict_event[f'semi_inclusive_chjet_star_trigger_pt'].append(hadron.pt())
                    
                        # Search for recoil jets
                        for jet in jets_selected:
                            if abs(jet.eta()) < (self.semi_inclusive_chjet_observables['IAA_star']['eta_cut_R'] - jetR):
                                            
                                # Get the corrected jet pt: shower+recoil-holes
                                jet_pt = jet.pt()
                                for temp_hadron in fj_hadrons_negative_charged:
                                    if jet.delta_R(temp_hadron) < jetR:
                                        jet_pt -= temp_hadron.pt()

                                # Jet yield and Delta phi
                                if jetR in self.semi_inclusive_chjet_observables['IAA_star']['jet_R']:

                                        if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                            if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_IAA_star_R{jetR}'].append(jet_pt)

                                        if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                            self.observable_dict_event[f'semi_inclusive_chjet_dphi_star_R{jetR}'].append(np.abs(hadron.delta_phi_to(jet)))

    # ---------------------------------------------------------------
    # Fill dijet histograms
    # ---------------------------------------------------------------
    def fill_dijet_histograms(self, jets_selected_charged, fj_hadrons_negative_charged, jetR):

        if self.sqrts == 2760:
            x = 1

    # ---------------------------------------------------------------
    # Compute electric charge from pid
    # ---------------------------------------------------------------
    def charge(self, pid):
    
        if pid in [11, 13, -211, -321, -2212, -3222, 3112, 3312, 3334]:
            return -1.
        elif pid in [-11, -13, 211, 321, 2212, 3222, -3112, -3312, -3334]:
            return 1.
        elif pid in [22, 111, 2112]:
            return 0.
        else:
            sys.exit(f'failed to compute charge of pid {pid}')

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
