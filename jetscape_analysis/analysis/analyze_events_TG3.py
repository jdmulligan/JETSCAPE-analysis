#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE output file

  Author: James Mulligan (james.mulligan@berkeley.edu)
  Author: Raymond Ehlers (raymond.ehlers@cern.ch)
  """

from __future__ import print_function

# General
import sys
import os
import argparse
import yaml
import numpy as np
import random
from collections import defaultdict

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
        self.constituent_threshold = config['constituent_threshold']
    
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
        self.inclusive_chjet_angularity_alice_bins = defaultdict(dict)
        for label in ['groomed', 'ungroomed']:
            self.inclusive_chjet_angularity_alice_bins[label] = defaultdict(dict)
            for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                for i,_ in enumerate(self.inclusive_chjet_observables['angularity_alice']['pt']):
                    if i < len(self.inclusive_chjet_observables['angularity_alice']['pt'])-1:
                        pt_min = self.inclusive_chjet_observables['angularity_alice']['pt'][i]
                        pt_max = self.inclusive_chjet_observables['angularity_alice']['pt'][i+1]
                        self.inclusive_chjet_angularity_alice_bins[label][f'{pt_min}-{pt_max}'][alpha] = np.array(self.inclusive_chjet_observables['angularity_alice']['bins'][label][f'{pt_min}-{pt_max}'][alpha])

        # Jet mass (binning not in hepdata)
        self.inclusive_chjet_mass_alice_bins = np.array(self.inclusive_chjet_observables['mass_alice']['bins'])

        # Soft Drop
        self.inclusive_chjet_zg_alice_bins = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_zg_central'])
        self.inclusive_chjet_tg_alice_bins = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_tg_central'])
        self.inclusive_chjet_zg_alice_bins_semicentral = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_zg_semicentral'])
        self.inclusive_chjet_tg_alice_bins_semicentral_zcut02 = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_tg_semicentral_zcut02'])
        self.inclusive_chjet_tg_alice_bins_semicentral_zcut04 = np.array(self.inclusive_chjet_observables['softdrop_alice']['bins_tg_semicentral_zcut04'])

        # Subjet z
        self.inclusive_chjet_subjets_alice_bins = np.linspace(0., 1., 100+1)
        
        # Jet axis
        self.inclusive_chjet_axis_alice_bins = {}
        self.inclusive_chjet_axis_alice_bins['Standard_WTA_R0.2'] = np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_WTA'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut01_R0.2'] = np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut02_R0.2'] = np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut02'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut01_R0.2'] = np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut02_R0.2'] = np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut02'])

        self.inclusive_chjet_axis_alice_bins['Standard_WTA_R0.4'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_WTA'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut01_R0.4'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut02_R0.4'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut02'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut01_R0.4'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut02_R0.4'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut02'])

        self.inclusive_chjet_axis_alice_bins['Standard_WTA_R0.5'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_WTA'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut01_R0.5'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['WTA_SD_zcut02_R0.5'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_WTA_SD_zcut02'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut01_R0.5'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut01'])
        self.inclusive_chjet_axis_alice_bins['Standard_SD_zcut02_R0.5'] = 2*np.array(self.inclusive_chjet_observables['axis_alice']['bins_Standard_SD_zcut02'])

        # Hardest kt
        self.inclusive_chjet_hardest_kt_alice_R02_bins = np.array(self.inclusive_chjet_observables["hardest_kt_alice"]["bins_ktg_R02"])
        self.inclusive_chjet_hardest_kt_alice_R04_bins = np.array(self.inclusive_chjet_observables["hardest_kt_alice"]["bins_ktg_R04"])
        self.inclusive_chjet_hardest_kt_alice_R05_bins = np.array(self.inclusive_chjet_observables["hardest_kt_alice"]["bins_ktg_R05"])
        
        # Charged jet RAA
        self.inclusive_chjet_pt_bins = np.linspace(0., 300., 300+1)
        
        #------------------------------------------------------
        # Semi-inclusive jet binnings
        
        # Hadron-jet IAA, dphi
        f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_IAA_276_R04'], 'READ')
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
        if self.user_index_for_pid:
            self.initialize_hadron_histograms()
        self.initialize_inclusive_jet_histograms()
        self.initialize_inclusive_chjet_histograms()
        self.initialize_semi_inclusive_chjet_histograms(charged=True)
        self.initialize_semi_inclusive_chjet_histograms(charged=False)

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
        
        # Holes
        hname = 'h_hadron_pt_cms_holes'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_cms_bins)-1, self.hadron_pt_cms_bins)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'h_hadron_pt_atlas_holes'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_atlas_bins)-1, self.hadron_pt_atlas_bins)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'h_hadron_pt_alice_holes'
        h = ROOT.TH1F(hname, hname, len(self.hadron_pt_alice_bins)-1, self.hadron_pt_alice_bins)
        h.Sumw2()
        setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_inclusive_jet_histograms(self):
    
        # Inclusive full jet histograms
        for jetR in self.jet_R:
            for constituent_threshold in self.constituent_threshold:

        
                hname = f'h_jet_pt_atlas_R{jetR}_0_10_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_0_10)-1,
                                                self.inclusive_jet_pt_atlas_bins_0_10)
                h.Sumw2()
                setattr(self, hname, h)
                
                hname = f'h_jet_pt_atlas_R{jetR}_30_40_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_30_40)-1,
                                                self.inclusive_jet_pt_atlas_bins_30_40)
                h.Sumw2()
                setattr(self, hname, h)
                
                hname = f'h_jet_pt_atlas_R{jetR}_40_50_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_atlas_bins_40_50)-1,
                                                self.inclusive_jet_pt_atlas_bins_40_50)
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_jet_pt_cms_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_cms_bins)-1,
                                                self.inclusive_jet_pt_cms_bins)
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_jet_pt_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_alice)-1,
                                                self.inclusive_jet_pt_alice)
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_jet_pt_alice_no_ptlead_cut_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_jet_pt_alice)-1,
                                            self.inclusive_jet_pt_alice)
                h.Sumw2()
                setattr(self, hname, h)
                
                hname = f'h_jet_pt_recoils_R{jetR}_pt{constituent_threshold}'
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
            for constituent_threshold in self.constituent_threshold:

                # g
                hname = f'h_chjet_g_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_g_alice_bins)-1,
                                                self.inclusive_chjet_g_alice_bins)
                h.Sumw2()
                setattr(self, hname, h)
                
                # Angularity (5.02 definition)
                for label in ['groomed', 'ungroomed']:
                    for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                        for i,_ in enumerate(self.inclusive_chjet_observables['angularity_alice']['pt']):
                            if i < len(self.inclusive_chjet_observables['angularity_alice']['pt'])-1:
                                pt_min = self.inclusive_chjet_observables['angularity_alice']['pt'][i]
                                pt_max = self.inclusive_chjet_observables['angularity_alice']['pt'][i+1]
                                hname = f'h_chjet_angularity_{label}_alice_R{jetR}_pt{pt_min}-{pt_max}_alpha{alpha}_pt{constituent_threshold}'
                                h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_angularity_alice_bins[label][f'{pt_min}-{pt_max}'][alpha])-1, 
                                                            self.inclusive_chjet_angularity_alice_bins[label][f'{pt_min}-{pt_max}'][alpha])
                                h.Sumw2()
                                setattr(self, hname, h)

                # Jet mass
                hname = f'h_chjet_mass_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_mass_alice_bins)-1,
                                                self.inclusive_chjet_mass_alice_bins)
                h.Sumw2()
                setattr(self, hname, h)

                # Soft Drop
                for zcut in self.inclusive_chjet_observables['soft_drop_zcut']:
                    if np.isclose(jetR, 0.2):
                        zg_bins_central = self.inclusive_chjet_zg_alice_bins
                        zg_bins_semicentral = self.inclusive_chjet_zg_alice_bins
                        tg_bins = self.inclusive_chjet_tg_alice_bins
                    else:
                        zg_bins_central = self.inclusive_chjet_zg_alice_bins
                        zg_bins_semicentral = self.inclusive_chjet_zg_alice_bins_semicentral
                        if np.isclose(zcut, 0.2):
                            tg_bins = self.inclusive_chjet_tg_alice_bins_semicentral_zcut02
                        else:
                            tg_bins = self.inclusive_chjet_tg_alice_bins_semicentral_zcut04
                    
                    hname = f'h_chjet_zg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}_central'
                    h = ROOT.TH1F(hname, hname, len(zg_bins_central)-1, zg_bins_central)
                    h.Sumw2()
                    setattr(self, hname, h)
                    
                    hname = f'h_chjet_zg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}_semicentral'
                    h = ROOT.TH1F(hname, hname, len(zg_bins_semicentral)-1, zg_bins_semicentral)
                    h.Sumw2()
                    setattr(self, hname, h)

                    hname = f'h_chjet_tg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}'
                    h = ROOT.TH1F(hname, hname, len(tg_bins)-1, tg_bins)
                    h.Sumw2()
                    setattr(self, hname, h)
                
                # Subjet z
                for r in self.inclusive_chjet_observables['subjetz_alice']['r']:
                    hname = f'h_chjet_subjetz_alice_R{jetR}_r{r}_pt{constituent_threshold}'
                    h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins,
                                                len(self.inclusive_chjet_subjets_alice_bins)-1,
                                                self.inclusive_chjet_subjets_alice_bins)
                    h.Sumw2()
                    setattr(self, hname, h)
            
                # Jet axis
                hname = f'h_chjet_axis_Standard_WTA_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins, 
                                                len(self.inclusive_chjet_axis_alice_bins[f'Standard_WTA_R{jetR}'])-1,
                                                self.inclusive_chjet_axis_alice_bins[f'Standard_WTA_R{jetR}'])
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_chjet_axis_WTA_SD_alice_R{jetR}_zcut01_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins, 
                                                len(self.inclusive_chjet_axis_alice_bins[f'WTA_SD_zcut01_R{jetR}'])-1,
                                                self.inclusive_chjet_axis_alice_bins[f'WTA_SD_zcut01_R{jetR}'])
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_chjet_axis_WTA_SD_alice_R{jetR}_zcut02_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins, 
                                                len(self.inclusive_chjet_axis_alice_bins[f'WTA_SD_zcut02_R{jetR}'])-1,
                                                self.inclusive_chjet_axis_alice_bins[f'WTA_SD_zcut02_R{jetR}'])
                h.Sumw2()
                setattr(self, hname, h)
                
                hname = f'h_chjet_axis_Standard_SD_alice_R{jetR}_zcut01_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins, 
                                                len(self.inclusive_chjet_axis_alice_bins[f'Standard_SD_zcut01_R{jetR}'])-1,
                                                self.inclusive_chjet_axis_alice_bins[f'Standard_SD_zcut01_R{jetR}'])
                h.Sumw2()
                setattr(self, hname, h)
                
                hname = f'h_chjet_axis_Standard_SD_alice_R{jetR}_zcut02_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, len(self.inclusive_chjet_pt_bins)-1, self.inclusive_chjet_pt_bins, 
                                                len(self.inclusive_chjet_axis_alice_bins[f'Standard_SD_zcut02_R{jetR}'])-1,
                                                self.inclusive_chjet_axis_alice_bins[f'Standard_SD_zcut02_R{jetR}'])
                h.Sumw2()
                setattr(self, hname, h)

                # Hardest kt
                for a in self.inclusive_chjet_observables["hardest_kt_alice"]["dynamical_grooming_a"]:
                    hname = f"h_chjet_ktg_dyg_a_{round(a*10):03}_alice_R{round(jetR*10):02}_pt{constituent_threshold}"
                    h = ROOT.TH1F(hname, hname, len(getattr(self, f"inclusive_chjet_hardest_kt_alice_R{round(jetR * 10):02}_bins"))-1,
                                                getattr(self, f"inclusive_chjet_hardest_kt_alice_R{round(jetR * 10):02}_bins"))
                    h.Sumw2()
                    setattr(self, hname, h)

                hname = f"h_chjet_ktg_soft_drop_z_cut_02_alice_R{round(jetR*10):02}_pt{constituent_threshold}"
                h = ROOT.TH1F(hname, hname, len(getattr(self, f"inclusive_chjet_hardest_kt_alice_R{round(jetR * 10):02}_bins"))-1,
                                            getattr(self, f"inclusive_chjet_hardest_kt_alice_R{round(jetR * 10):02}_bins"))
                h.Sumw2()
                setattr(self, hname, h)

                hname = f'h_chjet_pt_recoils_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH2F(hname, hname, 100, 0, 1000, 1000, 0, 100)
                h.GetXaxis().SetTitle('chjet pt')
                h.GetYaxis().SetTitle('recoil pt')
                h.Sumw2()
                setattr(self, hname, h)
                
                # Charged jet pt
                hname = f'h_chjet_pt_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(self.inclusive_chjet_pt_bins)-1,
                                                self.inclusive_chjet_pt_bins)
                h.Sumw2()
                setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_semi_inclusive_chjet_histograms(self, charged=True):

        if charged:
            ch_label = 'chjet'
        else:
            ch_label = 'jet'

        for jetR in self.jet_R:
            for constituent_threshold in self.constituent_threshold:

                # h-jet
                for hist_label in ['low', 'high']:
                
                    # Yield
                    hname = f'h_semi_inclusive_{ch_label}_IAA_{hist_label}Trigger_alice_R{jetR}_276_pt{constituent_threshold}'
                    h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_IAA_alice_276_bins)-1,
                                                    self.semi_inclusive_chjet_IAA_alice_276_bins)
                    h.Sumw2()
                    setattr(self, hname, h)
                    
                    # Delta Phi
                    hname = f'h_semi_inclusive_{ch_label}_dphi_{hist_label}Trigger_alice_R{jetR}_276_pt{constituent_threshold}'
                    h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_dphi_alice_276_bins)-1,
                                                    self.semi_inclusive_chjet_dphi_alice_276_bins)
                    h.Sumw2()
                    setattr(self, hname, h)

                    # For 2.76 TeV, make the 2D hist for folding (but keep the above so we can make predictions easily)
                    hname = f'h_semi_inclusive_{ch_label}_IAA_dphi_{hist_label}Trigger_alice_R{jetR}_276_pt{constituent_threshold}'
                    h = ROOT.TH2F(hname, hname, len(self.semi_inclusive_chjet_IAA_alice_502_bins)-1,
                                                self.semi_inclusive_chjet_IAA_alice_502_bins,
                                                len(self.semi_inclusive_chjet_dphi_alice_502_bins)-1,
                                                self.semi_inclusive_chjet_dphi_alice_502_bins)
                    h.Sumw2()
                    setattr(self, hname, h)
                    
                    # For 5.02 TeV, make 2D hist instead
                    hname = f'h_semi_inclusive_{ch_label}_IAA_dphi_{hist_label}Trigger_alice_R{jetR}_502_pt{constituent_threshold}'
                    h = ROOT.TH2F(hname, hname, len(self.semi_inclusive_chjet_IAA_alice_502_bins)-1,
                                                self.semi_inclusive_chjet_IAA_alice_502_bins,
                                                len(self.semi_inclusive_chjet_dphi_alice_502_bins)-1,
                                                self.semi_inclusive_chjet_dphi_alice_502_bins)
                    h.Sumw2()
                    setattr(self, hname, h)
                    
                    # Nsubjettiness
                    hname = f'h_semi_inclusive_{ch_label}_nsubjettiness_{hist_label}Trigger_alice_R{jetR}_pt{constituent_threshold}'
                    h = ROOT.TH1F(hname, hname, len(self.semi_inclusive_chjet_nsubjettiness_alice_bins)-1,
                                                    self.semi_inclusive_chjet_nsubjettiness_alice_bins)
                    h.Sumw2()
                    setattr(self, hname, h)

                # N triggers
                bins = np.array([5., 7, 8, 9, 20, 50])
                hname = f'h_semi_inclusive_{ch_label}_hjet_ntrigger_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
                h.Sumw2()
                setattr(self, hname, h)
                
                # N triggers
                bins = np.array([8., 9, 15, 45])
                hname = f'h_semi_inclusive_{ch_label}_nsubjettiness_ntrigger_alice_R{jetR}_pt{constituent_threshold}'
                h = ROOT.TH1F(hname, hname, len(bins)-1, bins)
                h.Sumw2()
                setattr(self, hname, h)
                
    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    #
    # The jet finding is done on positive status particles (shower+recoil),
    # and the negative status particles (holes) using the negative recombiner
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Create list of fastjet::PseudoJets (separately for jet shower particles and holes)
        fj_hadrons_positive_all = self.fill_fastjet_constituents(event, select_status='+')
        fj_hadrons_negative_all = self.fill_fastjet_constituents(event, select_status='-')
        
        # Create list of charged particles
        fj_hadrons_positive_charged_all = self.fill_fastjet_constituents(event, select_status='+',
                                                                         select_charged=True)
        fj_hadrons_negative_charged_all = self.fill_fastjet_constituents(event, select_status='-',
                                                                         select_charged=True)

        # TODO: To use negative recombiner as well as pid, should make two sets of particles -- one where
        # user_index is set to status, and the other where it is set pid
        
        # Fill hadron histograms for jet shower particles
        if self.user_index_for_pid:
            self.fill_hadron_histograms(fj_hadrons_positive_all, status='+')
            self.fill_hadron_histograms(fj_hadrons_negative_all, status='-')

        # Loop through several different constituent thresholds
        for constituent_threshold in self.constituent_threshold:
                    
            # Set constituent threshold
            fj_hadrons_positive = [hadron for hadron in fj_hadrons_positive_all if hadron.pt() > constituent_threshold]
            fj_hadrons_negative = [hadron for hadron in fj_hadrons_negative_all if hadron.pt() > constituent_threshold]
            fj_hadrons_positive_charged = [hadron for hadron in fj_hadrons_positive_charged_all if hadron.pt() > constituent_threshold]
            fj_hadrons_negative_charged = [hadron for hadron in fj_hadrons_negative_charged_all if hadron.pt() > constituent_threshold]

            if self.is_AA:
                hadrons_for_jet_finding = list(fj_hadrons_positive) + list(fj_hadrons_negative)
                hadrons_for_jet_finding_charged = list(fj_hadrons_positive_charged) + list(fj_hadrons_negative_charged)
            else:
                hadrons_for_jet_finding = list(fj_hadrons_positive)
                hadrons_for_jet_finding_charged = list(fj_hadrons_positive_charged)
                
            # Loop through specified jet R
            for jetR in self.jet_R:
                
                # Set jet definition and a jet selector
                jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
                if self.is_AA:
                    recombiner = fjext.NegativeEnergyRecombiner()
                    jet_def.set_recombiner(recombiner)
                jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(self.max_jet_y)
                if self.debug_level > 0:
                    print('jet definition is:', jet_def)
                    print('jet selector is:', jet_selector, '\n')

                # Full jets
                # -----------------
                cs = fj.ClusterSequence(hadrons_for_jet_finding, jet_def)
                jets = fj.sorted_by_pt(cs.inclusive_jets())
                jets_selected = jet_selector(jets)

                # Fill inclusive full jet histograms
                [self.analyze_inclusive_jet(jet, fj_hadrons_positive, fj_hadrons_negative, jetR, constituent_threshold, charged=False) for jet in jets_selected]
                
                # Fill jet correlations -- full
                self.fill_semi_inclusive_jet_histograms(jets_selected, fj_hadrons_positive, fj_hadrons_negative, jetR, constituent_threshold, charged=False)
                
                # Charged jets
                # -----------------
                cs_charged = fj.ClusterSequence(hadrons_for_jet_finding_charged, jet_def)
                jets_charged = fj.sorted_by_pt(cs_charged.inclusive_jets())
                jets_selected_charged = jet_selector(jets_charged)

                # Fill inclusive charged jet histograms
                [self.analyze_inclusive_jet(jet, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR, constituent_threshold, charged=True) for jet in jets_selected_charged]
                
                # Fill jet correlations -- charged
                self.fill_semi_inclusive_jet_histograms(jets_selected_charged, fj_hadrons_positive_charged, fj_hadrons_negative_charged, jetR, constituent_threshold, charged=True)

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

            # CMS
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.hadron_observables['pt_cms']['eta_cut']:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    if status == '-':
                        getattr(self, 'h_hadron_pt_cms_holes').Fill(pt)
                    else:
                        getattr(self, 'h_hadron_pt_cms').Fill(pt)

            # ATLAS
            # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            # Exclude e+, mu+ (11, 13)
            if abs(eta) < self.hadron_observables['pt_atlas']['eta_cut']:
                if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                    if status == '-':
                        getattr(self, 'h_hadron_pt_atlas_holes').Fill(pt)
                    else:
                        getattr(self, 'h_hadron_pt_atlas').Fill(pt)
                        
            # ALICE
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.hadron_observables['pt_alice']['eta_cut']:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    if status == '-':
                        getattr(self, 'h_hadron_pt_alice_holes').Fill(pt)
                    else:
                        getattr(self, 'h_hadron_pt_alice').Fill(pt)

    # ---------------------------------------------------------------
    # Fill inclusive jet histograms
    #
    # To correct jet pt: taken care of by negative recombiner
    # To correct substructure:
    #   - For additive observables, sum up the hole substructure observable and subtract
    #   - For identified objects within jets or between jets (groomed jet, subjet, jet axis, delta_phi),
    #     construct the observable only from the negative recombiner
    # ---------------------------------------------------------------
    def analyze_inclusive_jet(self, jet, fj_hadrons_positive, fj_hadrons_negative, jetR, constituent_threshold, charged=False):

        # Get the list of holes inside the jet -- we do not need to adjust the pt, but we want to keep track of the holes
        holes_in_jet = []
        negative_pt = 0.
        for hadron in fj_hadrons_negative:
            if jet.delta_R(hadron) < jetR:
                negative_pt += hadron.pt()
                holes_in_jet.append(hadron)
    
        jet_pt = jet.pt()
        
        if charged:
            getattr(self, f'h_chjet_pt_recoils_R{jetR}_pt{constituent_threshold}').Fill(jet_pt, negative_pt)
        else:
            getattr(self, f'h_jet_pt_recoils_R{jetR}_pt{constituent_threshold}').Fill(jet_pt, negative_pt)
            
        # Construct groomed jet
        # For negative_recombiner case, we set the negative recombiner also for the C/A reclustering
        jet_def = fj.JetDefinition(fj.cambridge_algorithm, jetR)
        if self.is_AA:
            recombiner = fjext.NegativeEnergyRecombiner()
            jet_def.set_recombiner(recombiner)
        gshop = fjcontrib.GroomerShop(jet, jet_def)

        # Fill histograms
        if charged:
            self.fill_charged_jet_histograms(jet, gshop, holes_in_jet, jet_pt, jetR, constituent_threshold)
        else:
            self.fill_full_jet_histograms(jet, jet_pt, jetR, constituent_threshold)

    # ---------------------------------------------------------------
    # Fill inclusive full jet histograms
    # ---------------------------------------------------------------
    def fill_full_jet_histograms(self, jet, jet_pt, jetR, constituent_threshold):

        # CMS RAA
        if abs(jet.eta()) < self.inclusive_jet_observables['pt_cms']['eta_cut']:
            getattr(self, f'h_jet_pt_cms_R{jetR}_pt{constituent_threshold}').Fill(jet_pt)

        # ATLAS RAA
        if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
            getattr(self, f'h_jet_pt_atlas_R{jetR}_0_10_pt{constituent_threshold}').Fill(jet_pt)
            getattr(self, f'h_jet_pt_atlas_R{jetR}_30_40_pt{constituent_threshold}').Fill(jet_pt)
            getattr(self, f'h_jet_pt_atlas_R{jetR}_40_50_pt{constituent_threshold}').Fill(jet_pt)

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
                getattr(self, f'h_jet_pt_alice_R{jetR}_pt{constituent_threshold}').Fill(jet_pt)
            getattr(self, f'h_jet_pt_alice_no_ptlead_cut_R{jetR}_pt{constituent_threshold}').Fill(jet_pt)
                
    # ---------------------------------------------------------------
    # Fill inclusive charged jet histograms
    # ---------------------------------------------------------------
    def fill_charged_jet_histograms(self, jet, gshop, holes_in_jet, jet_pt, jetR, constituent_threshold):
    
        # Grooming jet
        jet_groomed_lund = {0.1: None, 0.2: None}
        for zcut in self.inclusive_chjet_observables['soft_drop_zcut']:
            jet_groomed_lund[zcut] = gshop.soft_drop(self.inclusive_chjet_observables['soft_drop_beta'], zcut, jetR)
    
        # g
        if 40 < jet_pt < 60 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            g = 0
            for hadron in jet.constituents():
                if hadron.user_index() > 0:
                    g += hadron.pt() / jet_pt * hadron.delta_R(jet)
                for hadron in holes_in_jet:
                    if hadron.user_index() > 0:
                        continue
                    g -= hadron.pt() / jet_pt * hadron.delta_R(jet)
            getattr(self, f'h_chjet_g_alice_R{jetR}_pt{constituent_threshold}').Fill(g)

        # Angularity (5.02 definition)
        if 20 < jet_pt < 150 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                for i,_ in enumerate(self.inclusive_chjet_observables['angularity_alice']['pt']):
                    if i < len(self.inclusive_chjet_observables['angularity_alice']['pt'])-1:
                        pt_min = self.inclusive_chjet_observables['angularity_alice']['pt'][i]
                        pt_max = self.inclusive_chjet_observables['angularity_alice']['pt'][i+1]
                        if pt_min < jet_pt < pt_max:

                            lambda_alpha = 0
                            for hadron in jet.constituents():
                                if hadron.user_index() > 0:
                                    lambda_alpha += hadron.pt() / jet_pt * np.power(hadron.delta_R(jet)/jetR, alpha)
                                elif hadron.user_index() < 0:
                                    lambda_alpha -= hadron.pt() / jet_pt * np.power(hadron.delta_R(jet)/jetR, alpha)
                            getattr(self, f'h_chjet_angularity_ungroomed_alice_R{jetR}_pt{pt_min}-{pt_max}_alpha{alpha}_pt{constituent_threshold}').Fill(lambda_alpha)
                            
                            if jet_groomed_lund[0.2]:
                                lambda_alpha_g = 0
                                if jet_groomed_lund[0.2].pair().has_constituents():
                                    for hadron in jet_groomed_lund[0.2].pair().constituents():
                                        if hadron.user_index() > 0:
                                            lambda_alpha_g += hadron.pt() / jet_pt * np.power(hadron.delta_R(jet)/jetR, alpha)
                                        elif hadron.user_index() < 0:
                                            lambda_alpha_g -= hadron.pt() / jet_pt * np.power(hadron.delta_R(jet)/jetR, alpha)
                                    getattr(self, f'h_chjet_angularity_groomed_alice_R{jetR}_pt{pt_min}-{pt_max}_alpha{alpha}_pt{constituent_threshold}').Fill(lambda_alpha_g)

        # Jet mass
        if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            jet_mass = jet.m()
            getattr(self, f'h_chjet_mass_alice_R{jetR}_pt{constituent_threshold}').Fill(jet_mass)
            
        # Subjet z
        if 20 < jet_pt < 200 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            for r in self.inclusive_chjet_observables['subjetz_alice']['r']:
                
                cs_subjet = fj.ClusterSequence(jet.constituents(), fj.JetDefinition(fj.antikt_algorithm, r))
                subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())

                # Get leading subjet (accounts for holes)
                leading_subjet, leading_subjet_pt = self.leading_jet(subjets, holes_in_jet, r)
                z_leading = leading_subjet_pt / jet_pt
                
                # If z=1, it will be default be placed in overflow bin -- prevent this
                if np.isclose(z_leading, 1.):
                    z_leading = 0.999
                
                # Fill histogram
                getattr(self, f'h_chjet_subjetz_alice_R{jetR}_r{r}_pt{constituent_threshold}').Fill(jet_pt, z_leading)
        
        # Jet axis
        if 20 < jet_pt < 300 and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            
            # Recluster with WTA (with larger jet R)
            jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
            if self.is_AA:
                recombiner = fjext.NegativeEnergyRecombiner()
                jet_def_wta.set_recombiner(recombiner)
            jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
            reclusterer_wta = fjcontrib.Recluster(jet_def_wta)
            jet_wta = reclusterer_wta.result(jet)
            
            # Standard-WTA
            deltaR = jet.delta_R(jet_wta)
            getattr(self, f'h_chjet_axis_Standard_WTA_alice_R{jetR}_pt{constituent_threshold}').Fill(jet_pt, deltaR)

            # SD-WTA
            if jet_groomed_lund[0.1]:
                deltaR = jet_wta.delta_R(jet_groomed_lund[0.1].pair())
                getattr(self, f'h_chjet_axis_WTA_SD_alice_R{jetR}_zcut01_pt{constituent_threshold}').Fill(jet_pt, deltaR)
            if jet_groomed_lund[0.2]:
                deltaR = jet_wta.delta_R(jet_groomed_lund[0.2].pair())
                getattr(self, f'h_chjet_axis_WTA_SD_alice_R{jetR}_zcut02_pt{constituent_threshold}').Fill(jet_pt, deltaR)

            # Standard-SD
            if jet_groomed_lund[0.1]:
                deltaR = jet.delta_R(jet_groomed_lund[0.1].pair())
                getattr(self, f'h_chjet_axis_Standard_SD_alice_R{jetR}_zcut01_pt{constituent_threshold}').Fill(jet_pt, deltaR)
            if jet_groomed_lund[0.2]:
                deltaR = jet.delta_R(jet_groomed_lund[0.2].pair())
                getattr(self, f'h_chjet_axis_Standard_SD_alice_R{jetR}_zcut02_pt{constituent_threshold}').Fill(jet_pt, deltaR)

        # Hardest kt
        if 60 < jet_pt < 80 and abs(jet.eta()) < (self.inclusive_chjet_observables["eta_cut_alice_R"] - jetR):
            for a in self.inclusive_chjet_observables["hardest_kt_alice"]["dynamical_grooming_a"]:
                groomed = gshop.dynamical(a)
                ktg = groomed.kt()
                # Note: untagged will return kt = 0
                getattr(self, f"h_chjet_ktg_dyg_a_{round(a*10):03}_alice_R{round(jetR*10):02}_pt{constituent_threshold}").Fill(ktg)

            ktg = jet_groomed_lund[0.2].kt()
            # Note: untagged jets will return kt = 0
            getattr(self, f"h_chjet_ktg_soft_drop_z_cut_02_alice_R{round(jetR*10):02}_pt{constituent_threshold}").Fill(ktg)
            
        # Charged jet pt
        if abs(jet.eta()) < (self.inclusive_jet_observables['pt_alice']['eta_cut_R'] - jetR):
            getattr(self, f'h_chjet_pt_R{jetR}_pt{constituent_threshold}').Fill(jet_pt)


        for zcut in self.inclusive_chjet_observables['soft_drop_zcut']:
            self.fill_charged_groomed_jet_histograms(jet, gshop, jet_pt, jetR, constituent_threshold, zcut)

    # ---------------------------------------------------------------
    # Fill inclusive groomed charged jet histograms
    # ---------------------------------------------------------------
    def fill_charged_groomed_jet_histograms(self, jet, gshop, jet_pt, jetR, constituent_threshold, zcut):
    
        # Grooming jet
        jet_groomed_lund = gshop.soft_drop(self.inclusive_chjet_observables['soft_drop_beta'], zcut, jetR)
        
        if jet_groomed_lund and abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
            
            # Note: untagged jets will return negative value
            zg = jet_groomed_lund.z()
            theta_g = jet_groomed_lund.Delta() / jetR
        
            if 60 < jet_pt < 80:
                getattr(self, f'h_chjet_tg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}').Fill(theta_g)
        
            if np.isclose(jetR, 0.2):
                if 60 < jet_pt < 80:
                    getattr(self, f'h_chjet_zg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}_central').Fill(zg)
            else:
                if 60 < jet_pt < 80:
                    getattr(self, f'h_chjet_zg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}_semicentral').Fill(zg)
                if 80 < jet_pt < 100:
                    getattr(self, f'h_chjet_zg_alice_R{jetR}_pt{constituent_threshold}_zcut{zcut}_central').Fill(zg)

    # ---------------------------------------------------------------
    # Fill semi-inclusive charged jet histograms
    #
    # Note: We may need a lower jet pt range to determine cref, but I didn't look into this.
    # Note: Doesn't account for detector effects on hadron.
    # ---------------------------------------------------------------
    def fill_semi_inclusive_jet_histograms(self, jets_selected, fj_hadrons_positive_charged,
                                             fj_hadrons_negative_charged, jetR, constituent_threshold, charged=True):

        if charged:
            ch_label = 'chjet'
        else:
            ch_label = 'jet'

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

        # split events into signal and reference-classed events
        # majority signal to optimise stat. unc.
        frac_signal = 0.8
        is_signal_event = True
        if random.random() > frac_signal:
            is_signal_event = False


        # hadron + jet
        trigger_array_hjet = []

        for hadron in fj_hadrons_positive_charged:
        
            if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['hjet_alice']['hadron_eta_cut']:

                # Search for hadron trigger
                if hjet_low_trigger_range_276[0] < hadron.pt() < hjet_low_trigger_range_276[1] and not is_signal_event:
                    trigger_array_hjet.append(hadron)
                if hjet_high_trigger_range[0] < hadron.pt() < hjet_high_trigger_range[1] and is_signal_event:
                    trigger_array_hjet.append(hadron)

        # Search for recoil jets
        if len(trigger_array_hjet) > 0:

            # random selection of the trigger, since we may have more than one found in the event
            trigger = trigger_array_hjet[random.randrange(len(trigger_array_hjet))]

            # Record hadron pt for trigger normalization
            # NOTE: This will record the hadron trigger even if it's not used in the IAA. However,
            #       this is fine because we account for the difference in low and high trigger ranges
            #       when we construct the histograms.
            getattr(self, f'h_semi_inclusive_{ch_label}_hjet_ntrigger_alice_R{jetR}_pt{constituent_threshold}').Fill(trigger.pt())

            for jet in jets_selected:
                if abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
                                
                    # Get the corrected jet pt: shower+recoil-holes
                    jet_pt = jet.pt()
                    for temp_hadron in fj_hadrons_negative_charged:
                        if jet.delta_R(temp_hadron) < jetR:
                            jet_pt -= temp_hadron.pt()

                    # Jet yield and Delta phi
                    if is_signal_event:
                        if np.abs(jet.delta_phi_to(trigger)) > (np.pi - 0.6):
                            getattr(self, f'h_semi_inclusive_{ch_label}_IAA_highTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(jet_pt)

                        if 40 < jet_pt < 60:
                            getattr(self, f'h_semi_inclusive_{ch_label}_dphi_highTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(np.abs(trigger.delta_phi_to(jet)))

                        getattr(self, f'h_semi_inclusive_{ch_label}_IAA_dphi_highTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(jet_pt, np.abs(trigger.delta_phi_to(jet)))

                    else:

                        if np.abs(jet.delta_phi_to(trigger)) > (np.pi - 0.6):
                            getattr(self, f'h_semi_inclusive_{ch_label}_IAA_lowTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(jet_pt)

                        if 40 < jet_pt < 60:
                            getattr(self, f'h_semi_inclusive_{ch_label}_dphi_lowTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(np.abs(trigger.delta_phi_to(jet)))

                        getattr(self, f'h_semi_inclusive_{ch_label}_IAA_dphi_lowTrigger_alice_R{jetR}_276_pt{constituent_threshold}').Fill(jet_pt, np.abs(trigger.delta_phi_to(jet)))

        # Nsubjettiness
        # We use the jet pt including recoils here, since the Nsubjettiness is calculated
        # including recoils (but without hole subtraction).
        # Not ideal but not sure of an immediate better solution.
                           
        trigger_array_nsubjettiness = []

        for hadron in fj_hadrons_positive_charged:
        
            if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['hjet_alice']['hadron_eta_cut']:

                # Search for hadron trigger
                if nsubjettiness_low_trigger_range[0] < hadron.pt() < nsubjettiness_low_trigger_range[1] and not is_signal_event:
                    trigger_array_nsubjettiness.append(hadron)
                if nsubjettiness_high_trigger_range[0] < hadron.pt() < nsubjettiness_high_trigger_range[1] and is_signal_event:
                    trigger_array_nsubjettiness.append(hadron)



        # Search for recoil jets
        if len(trigger_array_nsubjettiness) > 0:

            # random selection of the trigger, since we may have more than one found in the event
            trigger = trigger_array_nsubjettiness[random.randrange(len(trigger_array_nsubjettiness))]

            # Record hadron pt for trigger normalization
            getattr(self, f'h_semi_inclusive_{ch_label}_nsubjettiness_ntrigger_alice_R{jetR}_pt{constituent_threshold}').Fill(trigger.pt())

            for jet in jets_selected:
                if abs(jet.eta()) < (self.inclusive_chjet_observables['eta_cut_alice_R'] - jetR):
                                
                    # Get the corrected jet pt: shower+recoil-holes
                    jet_pt = jet.pt()
                    for temp_hadron in fj_hadrons_negative_charged:
                        if jet.delta_R(temp_hadron) < jetR:
                            jet_pt -= temp_hadron.pt()

                    if is_signal_event:
                        if np.abs(jet.delta_phi_to(trigger)) > (np.pi - 0.6):
                            if 40 < jet_pt < 60:
                                tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                if tau1 > 1e-3:
                                    getattr(self, f'h_semi_inclusive_{ch_label}_nsubjettiness_highTrigger_alice_R{jetR}_pt{constituent_threshold}').Fill(tau2/tau1)
                    else:
                        if np.abs(jet.delta_phi_to(trigger)) > (np.pi - 0.6):
                            if 40 < jet_pt < 60:
                                tau1 = n_subjettiness_calculator1.result(jet)/jet.pt()
                                tau2 = n_subjettiness_calculator2.result(jet)/jet.pt()
                                if tau1 > 1e-3:
                                    getattr(self, f'h_semi_inclusive_{ch_label}_nsubjettiness_lowTrigger_alice_R{jetR}_pt{constituent_threshold}').Fill(tau2/tau1)

    #---------------------------------------------------------------
    # Return leading jet (or subjet)
    #---------------------------------------------------------------
    def leading_jet(self, jets, fj_hadrons_negative, jetR):

        leading_jet = None
        leading_jet_pt = 0.
        for i,jet in enumerate(jets):
                         
            # Get the corrected jet pt by subtracting the negative recoils within R
            jet_pt = jet.pt()
            for temp_hadron in fj_hadrons_negative:
                if jet.delta_R(temp_hadron) < jetR:
                    jet_pt -= temp_hadron.pt()
                    
            if not leading_jet:
                leading_jet = jet
                leading_jet_pt = jet_pt
            
            if jet_pt > leading_jet_pt:
                leading_jet = jet
                leading_jet_pt = jet_pt

        return leading_jet, leading_jet_pt

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
