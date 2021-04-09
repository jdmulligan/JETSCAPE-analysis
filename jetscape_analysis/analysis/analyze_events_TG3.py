#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE output file

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

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import ROOT

sys.path.append('.')
from jetscape_analysis.analysis import analyze_events_base_PHYS

################################################################
class AnalyzeJetscapeEvents_TG3(analyze_events_base_PHYS.AnalyzeJetscapeEvents_BasePHYS):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_TG3, self).__init__(config_file=config_file,
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
            
        self.hadron_observables = config['hadron']
        self.inclusive_jet_observables = config['inclusive_jet']
        self.inclusive_chjet_observables = config['inclusive_chjet']
        self.semi_inclusive_chjet_observables = config['semi_inclusive_chjet']
        
        self.jet_R = config['jet_R']
        self.min_jet_pt = config['min_jet_pt']
        self.max_jet_y = config['max_jet_y']
    
    # ---------------------------------------------------------------
    # Initialize binnings
    # ---------------------------------------------------------------
    def initialize_binnings(self):

        #------------------------------------------------------
        # Charged particle binnings
        
        f = ROOT.TFile(self.hadron_observables['pt_cms']['hepdata_0_5'], 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.hadron_pt_cms_bins = np.array(h.GetXaxis().GetXbins())
        f.Close()

        self.hadron_pt_atlas_bins = np.array(self.hadron_observables['pt_atlas']['bins'])
        
        f = ROOT.TFile(self.hadron_observables['pt_alice']['hepdata'], 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.hadron_pt_alice_bins = np.array(h.GetXaxis().GetXbins())
        f.Close()

        #------------------------------------------------------
        # Inclusive full jet binnings

        self.inclusive_jet_pt_cms_bins = np.array(self.inclusive_jet_observables['pt_cms']['bins'])

        f = ROOT.TFile(self.inclusive_jet_observables['pt_atlas']['hepdata_0_10'], 'READ')
        dir = f.Get('Table 19')
        h = dir.Get('Hist1D_y1')
        self.inclusive_jet_pt_atlas_bins_0_10 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.inclusive_jet_observables['pt_atlas']['hepdata_30_40'], 'READ')
        dir = f.Get('Table 22')
        h = dir.Get('Hist1D_y1')
        self.inclusive_jet_pt_atlas_bins_30_40 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.inclusive_jet_observables['pt_atlas']['hepdata_40_50'], 'READ')
        dir = f.Get('Table 23')
        h = dir.Get('Hist1D_y1')
        self.inclusive_jet_pt_atlas_bins_40_50 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.inclusive_jet_observables['pt_alice']['hepdata_0_10_R02'], 'READ')
        dir = f.Get('Table 30')
        h = dir.Get('Hist1D_y1')
        self.inclusive_jet_pt_alice = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        #------------------------------------------------------
        # Inclusive charged jet binnings
        
        # g
        f = ROOT.TFile(self.inclusive_chjet_observables['g_alice']['hepdata'], 'READ')
        dir = f.Get('Table 11')
        h = dir.Get('Hist1D_y1')
        self.inclusive_chjet_g_alice_bins = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        # Angularity
        self.inclusive_chjet_angularity_alice_bins = np.linspace(0., 1., 100+1)

        # Jet mass (binning not in hepdata)
        self.inclusive_chjet_mass_alice_bins = np.array(self.inclusive_chjet_observables['mass_alice']['bins'])

        # Soft Drop
        self.inclusive_chjet_zg_alice_bins = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_zg'])
        self.inclusive_chjet_tg_alice_bins = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_tg'])
        
        # Subjet z
        self.inclusive_chjet_subjets_alice_bins = np.linspace(0., 1., 100+1)
        
        # Jet axis
        self.inclusive_chjet_axis_alice_bins = np.linspace(0., 0.2, 200+1)
        
        #------------------------------------------------------
        # Semi-inclusive jet binnings
        
        # Hadron-jet IAA, dphi
        f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_IAA_276'], 'READ')
        dir = f.Get('Table 33')
        h = dir.Get('Hist1D_y1')
        self.semi_inclusive_chjet_IAA_alice_276_bins = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_dphi_276'], 'READ')
        dir = f.Get('Table 37')
        h = dir.Get('Hist1D_y1')
        self.semi_inclusive_chjet_dphi_alice_276_bins = np.array(h.GetXaxis().GetXbins())
        f.Close()
        
        self.semi_inclusive_chjet_IAA_alice_502_bins = np.linspace(0., 200., 200+1)
        self.semi_inclusive_chjet_dphi_alice_502_bins = np.linspace(np.pi/2, np.pi, 100+1)
        
        # Nsubjettiness
        self.semi_inclusive_chjet_nsubjettiness_alice_bins = np.linspace(0., 1., 100+1)
        
    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):
    
        # Construct binnings
        self.initialize_binnings()
        
        # Initialize each set of histograms
        self.initialize_hadron_histograms()
        self.initialize_inclusive_jet_histograms()
        self.initialize_inclusive_chjet_histograms()
        self.initialize_semi_inclusive_chjet_histograms()

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_hadron_histograms(self):
    
        # Hadron histograms
        hname = 'h_hadron_pt_cms'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_cms_bins)-1, self.hadron_pt_cms_bins)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'h_hadron_pt_atlas'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_atlas_bins)-1, self.hadron_pt_atlas_bins)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'h_hadron_pt_alice'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_alice_bins)-1, self.hadron_pt_alice_bins)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'h_hadron_pt_recoils'
        h = ROOT.TH1F(hname, hname, 1000, 0, 100)
        h.Sumw2()
        setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_inclusive_jet_histograms(self):
    
        # Inclusive full jet histograms
        for jetR in self.jet_R:
        
            hname = f'h_jet_pt_atlas_R{jetR}_0_10'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_0_10)-1,
                                            self.inclusive_jet_pt_atlas_bins_0_10)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_jet_pt_atlas_R{jetR}_30_40'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_30_40)-1,
                                            self.inclusive_jet_pt_atlas_bins_30_40)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_jet_pt_atlas_R{jetR}_40_50'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_40_50)-1,
                                            self.inclusive_jet_pt_atlas_bins_40_50)
            h.Sumw2()
            setattr(self, hname, h)

            hname = f'h_jet_pt_cms_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_cms_bins)-1,
                                            self.inclusive_jet_pt_cms_bins)
            h.Sumw2()
            setattr(self, hname, h)

            hname = f'h_jet_pt_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_alice)-1,
                                            self.inclusive_jet_pt_alice)
            h.Sumw2()
            setattr(self, hname, h)

            hname = f'h_jet_pt_alice_no_ptlead_cut_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_alice)-1,
                                        self.inclusive_jet_pt_alice)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_jet_pt_recoils_R{jetR}'
            h = ROOT.TH2F(hname, hname, 100, 0, 1000, 1000, 0, 100)
            h.GetXaxis().SetTitle('jet pt')
            h.GetYaxis().SetTitle('recoil pt')
            h.Sumw2()
            setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_inclusive_chjet_histograms(self):

        for jetR in self.jet_R:

            # g
            hname = f'h_chjet_g_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_g_alice_bins)-1,
                                            self.inclusive_chjet_g_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            # Angularity (5.02 definition)
            for label in ['groomed', 'ungroomed']:
                for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                    hname = f'h_chjet_angularity_{label}_alice_R{jetR}_alpha{alpha}'
                    h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_angularity_alice_bins)-1,
                                                self.inclusive_chjet_angularity_alice_bins)
                    h.Sumw2()
                    setattr(self, hname, h)

            # Jet mass
            hname = f'h_chjet_mass_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_mass_alice_bins)-1,
                                            self.inclusive_chjet_mass_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)

            # Soft Drop
            hname = f'h_chjet_zg_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_zg_alice_bins)-1,
                                            self.inclusive_chjet_zg_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)

            hname = f'h_chjet_tg_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_tg_alice_bins)-1,
                                        self.inclusive_chjet_tg_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            # Subjet z
            for r in self.inclusive_chjet_observables['subjetz_alice']['r']:
                hname = f'h_chjet_subjetz_alice_R{jetR}_r{r}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_subjets_alice_bins)-1,
                                                self.inclusive_chjet_subjets_alice_bins)
                h.Sumw2()
                setattr(self, hname, h)
        
            # Jet axis
            hname = f'h_chjet_axis_Standard_WTA_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_axis_alice_bins)-1,
                                            self.inclusive_chjet_axis_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_chjet_axis_Standard_SD_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_axis_alice_bins)-1,
                                            self.inclusive_chjet_axis_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_chjet_axis_SD_WTA_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_axis_alice_bins)-1,
                                            self.inclusive_chjet_axis_alice_bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            hname = f'h_chjet_pt_recoils_R{jetR}'
            h = ROOT.TH2F(hname, hname, 100, 0, 1000, 1000, 0, 100)
            h.GetXaxis().SetTitle('chjet pt')
            h.GetYaxis().SetTitle('recoil pt')
            h.Sumw2()
            setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_semi_inclusive_chjet_histograms(self):

        for jetR in self.jet_R:

            # h-jet
            for hist_label in ['low', 'high']:
            
                # Yield
                hname = f'h_semi_inclusive_chjet_IAA_{hist_label}Trigger_alice_R{jetR}_276'
                h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_IAA_alice_276_bins)-1,
                                                self.semi_inclusive_chjet_IAA_alice_276_bins)
                h.Sumw2()
                setattr(self, hname, h)
                
                # Delta Phi
                hname = f'h_semi_inclusive_chjet_dphi_{hist_label}Trigger_alice_R{jetR}_276'
                h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_dphi_alice_276_bins)-1,
                                                self.semi_inclusive_chjet_dphi_alice_276_bins)
                h.Sumw2()
                setattr(self, hname, h)
                
                # For 5.02 TeV, make 2D hist instead
                hname = f'h_semi_inclusive_chjet_IAA_dphi_{hist_label}Trigger_alice_R{jetR}_502'
                h = ROOT.TH2F(hname, hname, len(self.semi_inclusive_chjet_IAA_alice_502_bins)-1,
                                            self.semi_inclusive_chjet_IAA_alice_502_bins,
                                            len(self.semi_inclusive_chjet_dphi_alice_502_bins)-1,
                                            self.semi_inclusive_chjet_dphi_alice_502_bins)
                h.Sumw2()
                setattr(self, hname, h)
                
                # Nsubjettiness
                hname = f'h_semi_inclusive_chjet_nsubjettiness_{hist_label}Trigger_alice_R{jetR}'
                h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_nsubjettiness_alice_bins)-1,
                                                self.semi_inclusive_chjet_nsubjettiness_alice_bins)
                h.Sumw2()
                setattr(self, hname, h)

            # N triggers
            bins = np.array([5., 7, 8, 9, 20, 50])
            hname = f'h_semi_inclusive_chjet_hjet_ntrigger_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
            h.Sumw2()
            setattr(self, hname, h)
            
            # N triggers
            bins = np.array([8., 9, 15, 45])
            hname = f'h_semi_inclusive_chjet_nsubjettiness_ntrigger_alice_R{jetR}'
            h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
            h.Sumw2()
            setattr(self, hname, h)
            
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
        self.fill_hadron_histograms(fj_hadrons_negative, status='-')

        # Loop through specified jet R
        for jetR in self.jet_R:
        
            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(self.max_jet_y)
            if self.debug_level > 0:
                print('jet definition is:', jet_def)
                print('jet selector is:', jet_selector, '\n')

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
            
            # Fill jet correlations
            self.fill_semi_inclusive_chjet_histograms(jets_selected_charged, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR)

    # ---------------------------------------------------------------
    # Fill hadron histograms
    # (assuming weak strange decays are off, but charm decays are on)
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, fj_particles, status='+'):

        # Loop through hadrons
        for i,particle in enumerate(fj_particles):

            # Fill some basic hadron info
            pid = particle.user_index()
            pt = particle.pt()
            eta = particle.eta()

            # Skip negative recoils (holes)
            if status == '-':
                getattr(self, 'h_hadron_pt_recoils').Fill(pt)
                continue

            # CMS
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.hadron_observables['pt_cms']['eta_cut']:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'h_hadron_pt_cms').Fill(pt)

            # ATLAS
            # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            # Exclude e+, mu+ (11, 13)
            if abs(eta) < self.hadron_observables['pt_atlas']['eta_cut']:
                if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'h_hadron_pt_atlas').Fill(pt)

            # ALICE
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.hadron_observables['pt_alice']['eta_cut']:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'h_hadron_pt_alice').Fill(pt)

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
        
        if charged:
            getattr(self, f'h_chjet_pt_recoils_R{jetR}').Fill(jet_pt, negative_pt)
        else:
            getattr(self, f'h_jet_pt_recoils_R{jetR}').Fill(jet_pt, negative_pt)
            
        # Construct groomed jet
        gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
        jet_groomed_lund = gshop.soft_drop(self.inclusive_chjet_observables['soft_drop_beta'], self.inclusive_chjet_observables['soft_drop_zcut'], jetR)

        # Fill histograms
        if charged:
            self.fill_charged_jet_histograms(jet, jet_groomed_lund, holes_in_jet, jet_pt, jetR)
        else:
            self.fill_full_jet_histograms(jet, jet_pt, jetR)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_full_jet_histograms(self, jet, jet_pt, jetR):

        # CMS RAA
        if abs(jet.eta()) < self.inclusive_jet_observables['pt_cms']['eta_cut']:
            getattr(self, f'h_jet_pt_cms_R{jetR}').Fill(jet_pt)

        # ATLAS RAA
        if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
            getattr(self, f'h_jet_pt_atlas_R{jetR}_0_10').Fill(jet_pt)
            getattr(self, f'h_jet_pt_atlas_R{jetR}_30_40').Fill(jet_pt)
            getattr(self, f'h_jet_pt_atlas_R{jetR}_40_50').Fill(jet_pt)

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
                getattr(self, f'h_jet_pt_alice_R{jetR}').Fill(jet_pt)
            getattr(self, f'h_jet_pt_alice_no_ptlead_cut_R{jetR}').Fill(jet_pt)
    
    # ---------------------------------------------------------------
    # Fill inclusive charged jet histograms
    # ---------------------------------------------------------------
    def fill_charged_jet_histograms(self, jet, jet_groomed_lund, holes_in_jet, jet_pt, jetR):
    
        # g
        if 40 < jet_pt < 60 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            g = 0
            for constituent in jet.constituents():
                g += constituent.pt() / jet_pt * constituent.delta_R(jet)
            for hadron in holes_in_jet:
                g -= hadron.pt() / jet_pt * hadron.delta_R(jet)
            getattr(self, f'h_chjet_g_alice_R{jetR}').Fill(g)

        # Angularity (5.02 definition)
        if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
            
                kappa=1
                lambda_alpha = fjext.lambda_beta_kappa(jet, alpha, kappa, jetR)
                for hadron in holes_in_jet:
                    lambda_alpha -= hadron.pt() / jet_pt * np.power(hadron.delta_R(jet), alpha)
                getattr(self, f'h_chjet_angularity_ungroomed_alice_R{jetR}_alpha{alpha}').Fill(lambda_alpha)
                
                if jet_groomed_lund:
                    lambda_alpha_g = fjext.lambda_beta_kappa(jet_groomed_lund.pair(), alpha, kappa, jetR)
                    for hadron in holes_in_jet:
                        lambda_alpha_g -= hadron.pt() / jet_pt * np.power(hadron.delta_R(jet), alpha)
                    getattr(self, f'h_chjet_angularity_groomed_alice_R{jetR}_alpha{alpha}').Fill(lambda_alpha_g)

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

            getattr(self, f'h_chjet_mass_alice_R{jetR}').Fill(jet_mass)

        # Soft Drop
        if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            if jet_groomed_lund:
                theta_g = jet_groomed_lund.Delta() / jetR
                zg = jet_groomed_lund.z()
                # Note: untagged jets will return negative value
            getattr(self, f'h_chjet_zg_alice_R{jetR}').Fill(zg)
            getattr(self, f'h_chjet_tg_alice_R{jetR}').Fill(theta_g)
            
        # Subjet z
        if 80 < jet_pt < 100 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            for r in self.inclusive_chjet_observables['subjetz_alice']['r']:
                
                cs_subjet = fj.ClusterSequence(jet.constituents(), fj.JetDefinition(fj.antikt_algorithm, r))
                subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
                # Note: May be better to subtract holes before deciding leading subjet
                leading_subjet = self.leading_jet(subjets)
                
                # Sum the negative recoils within subjetR
                negative_pt = 0.
                for hadron in holes_in_jet:
                    if leading_subjet.delta_R(hadron) < r:
                        negative_pt += hadron.pt()

                # Compute corrected subjet pt, and fill histograms
                subjet_pt = leading_subjet.pt() - negative_pt
                z_leading = subjet_pt / jet_pt
                getattr(self, f'h_chjet_subjetz_alice_R{jetR}_r{r}').Fill(z_leading)
        
        # Jet axis
        if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            
            # Recluster with WTA (with larger jet R)
            jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
            jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
            reclusterer_wta = fjcontrib.Recluster(jet_def_wta)
            jet_wta = reclusterer_wta.result(jet)
            
            # Standard-WTA
            deltaR = jet.delta_R(jet_wta)
            getattr(self, f'h_chjet_axis_Standard_WTA_alice_R{jetR}').Fill(deltaR)

            # Standard-SD
            if jet_groomed_lund:
                deltaR = jet.delta_R(jet_groomed_lund.pair())
                getattr(self, f'h_chjet_axis_Standard_SD_alice_R{jetR}').Fill(deltaR)

            # SD-WTA
            if jet_groomed_lund:
                deltaR = jet_wta.delta_R(jet_groomed_lund.pair())
                getattr(self, f'h_chjet_axis_SD_WTA_alice_R{jetR}').Fill(deltaR)

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
                getattr(self, f'h_semi_inclusive_chjet_hjet_ntrigger_alice_R{jetR}').Fill(hadron.pt())
                getattr(self, f'h_semi_inclusive_chjet_nsubjettiness_ntrigger_alice_R{jetR}').Fill(hadron.pt())
                
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
                                    getattr(self, f'h_semi_inclusive_chjet_IAA_lowTrigger_alice_R{jetR}_276').Fill(jet_pt)

                                if 40 < jet_pt < 60:
                                    getattr(self, f'h_semi_inclusive_chjet_dphi_lowTrigger_alice_R{jetR}_276').Fill(np.abs(hadron.delta_phi_to(jet)))

                            if hjet_found_low_502:
                                getattr(self, f'h_semi_inclusive_chjet_IAA_dphi_lowTrigger_alice_R{jetR}_502').Fill(jet_pt, np.abs(hadron.delta_phi_to(jet)))
                                    
                            if hjet_found_high:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    getattr(self, f'h_semi_inclusive_chjet_IAA_highTrigger_alice_R{jetR}_276').Fill(jet_pt)

                                if 40 < jet_pt < 60:
                                    getattr(self, f'h_semi_inclusive_chjet_dphi_highTrigger_alice_R{jetR}_276').Fill(np.abs(hadron.delta_phi_to(jet)))

                                getattr(self, f'h_semi_inclusive_chjet_IAA_dphi_highTrigger_alice_R{jetR}_502').Fill(jet_pt, np.abs(hadron.delta_phi_to(jet)))

                            # Nsubjettiness
                            if nsubjettiness_found_low:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    if 40 < jet_pt < 60:
                                        tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                        tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                        if tau1 > 1e-3:
                                            getattr(self, f'h_semi_inclusive_chjet_nsubjettiness_lowTrigger_alice_R{jetR}').Fill(tau2/tau1)

                            if nsubjettiness_found_high:
                                if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                    if 40 < jet_pt < 60:
                                        tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                        tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                        if tau1 > 1e-3:
                                            getattr(self, f'h_semi_inclusive_chjet_nsubjettiness_highTrigger_alice_R{jetR}').Fill(tau2/tau1)

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

    analysis = AnalyzeJetscapeEvents_TG3(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
