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

# Data analysis and plotting
import ROOT
import ctypes
import numpy as np

# Base class
sys.path.append('.')
from jetscape_analysis.base import common_base
from plot import plot_results_STAT_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotResults(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', **kwargs):
        super(PlotResults, self).__init__(**kwargs)
        self.output_dir = os.path.dirname(input_file)
               
        self.plot_utils = plot_results_STAT_utils.PlotUtils()
        self.plot_utils.setOptions()
        ROOT.gROOT.ForceStyle()

        self.input_file = ROOT.TFile(input_file, 'READ')

        self.data_color = ROOT.kGray+3
        self.data_marker = 21
        self.jetscape_color = ROOT.kViolet-8
        self.jetscape_marker = 20
        self.alpha = 0.7
        self.marker_size = 1.5
        self.line_width = 2
        self.line_style = 1

        # Read config file
        with open(config_file, 'r') as stream:
            self.config = yaml.safe_load(stream)
        self.sqrts = self.config['sqrt_s']
      
        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):
    
        self.plot_hadron_observables(observable_type='hadron')
        
        #self.plot_hadron_correlation_observables(observable_type='hadron_correlations')
        
        #self.plot_jet_observables(observable_type='inclusive_chjet')
        
        #if 'inclusive_jet' in self.config:
        #    self.plot_jet_observables(observable_type='inclusive_jet')
            
        #if 'semi_inclusive_chjet' in self.config:
        #    self.plot_semi_inclusive_chjet_observables(observable_type='semi_inclusive_chjet')
            
        #if 'dijet' in self.config:
        #    self.plot_jet_observables(observable_type='dijet')

    #-------------------------------------------------------------------------------------------
    # Plot hadron observables
    #-------------------------------------------------------------------------------------------
    def plot_hadron_observables(self, observable_type=''):
        print()
        print(f'Plot {observable_type} observables...')
                
        for observable, block in self.config[observable_type].items():
        
            # Initialize observable configuration
            self.init_observable(observable_type, observable, block)
                
            # Plot observable
            self.plot_distribution_and_ratio(observable_type, observable, logy=True)

    #-------------------------------------------------------------------------------------------
    # Initialize a single observable's config
    #-------------------------------------------------------------------------------------------
    def init_observable(self, observable_type, observable, block, self_normalize=False):
    
        # Initialize an empty dict containing relevant info
        self.observable_settings = {}
                    
        # Common settings
        self.xtitle = block['xtitle']
        if 'eta_cut' in block:
            self.eta_cut = block['eta_cut']
        if 'pt' in block:
            self.pt = block['pt']
        if 'jetR' in block:
            self.jetR = block['jetR']
        if 'eta_R' in block:
            self.eta_R = block['eta_R']
            self.eta_cut = np.round(self.eta_R - self.jetR, decimals=1)
            
        if 'ytitle_pp' in block:
            self.ytitle = block['ytitle_pp']
        if 'y_min_pp' in block:
            self.y_min = float(block['y_min_pp'])
            self.y_max = float(block['y_max_pp'])
        if 'y_ratio_min' in block:
            self.y_ratio_min = block['y_ratio_min']
            self.y_ratio_max = block['y_ratio_max']
        else:
            self.y_ratio_min = 0.
            self.y_ratio_max = 1.99
        
        # Initialize data
        if 'hepdata' in block:
            f = ROOT.TFile(block['hepdata'], 'READ')
            dir = f.Get(block['hepdata_dir'])
            self.observable_settings['data_distribution'] = dir.Get(block['hepdata_gname'])
            if self.observable_settings['data_distribution'].InheritsFrom(ROOT.TH1.Class()):
                self.observable_settings['data_distribution'].SetDirectory(0)
            
        # Initialize JETSCAPE
        self.hname = f'h_{observable_type}_{observable}'
        h_jetscape = self.input_file.Get(self.hname)
        h_jetscape.SetDirectory(0)
        self.observable_settings['jetscape_distribution'] = h_jetscape
        
        n_events = self.input_file.Get('h_n_events').GetBinContent(1)
        if observable_type == 'hadron':
            self.observable_settings['jetscape_distribution'].Scale(1./n_events)
            self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
            sigma_inel = 0.0676
            self.observable_settings['jetscape_distribution'].Scale(1./sigma_inel)
        
        ## For hadrons, impose a 1 GeV minimum, and subtract the recoil hadrons
        #if self.observable == 'hadron_raa':

        #    # Subtract holes
        #    f_jetscape_AA = ROOT.TFile(block['file_jetscape_AA'], 'READ')
        #    h_recoil_name = 'h_hadron_pt_alice_holesScaled'
        #    h_recoil = f_jetscape_AA.Get(h_recoil_name)
        #    h_recoil.SetDirectory(0)
        #    h_recoil_rebinned = h_recoil.Rebin(h_jetscape_AA_xbins.size-1, f'{h_recoil_name}_AA', h_jetscape_AA_xbins)
        #    h_jetscape_AA.Add(h_recoil_rebinned, -1)
        #    f_jetscape_AA.Close()
            
        # Normalization
        h_jetscape.Scale(1., 'width')
        if self_normalize:
            if observable in ['zg', 'theta_g']:
                min_bin = 0
            else:
                min_bin = 1
            h_jetscape_pp.Scale(1./h_jetscape_pp.Integral(min_bin, h_jetscape_pp.GetNbinsX()))
            
        # Form ratio of JETSCAPE to data
        if self.observable_settings['data_distribution'].InheritsFrom(ROOT.TH1.Class()):
            h_ratio = self.observable_settings['jetscape_distribution']
            h_ratio.Divide(self.observable_settings['data_distribution'])
            self.observable_settings['ratio'] = h_ratio
        elif self.observable_settings['data_distribution'].InheritsFrom(ROOT.TGraph.Class()):
            self.observable_settings['ratio'] = self.plot_utils.divide_tgraph(self.observable_settings['jetscape_distribution'],
                                                                              self.observable_settings['data_distribution'])
                                                                              
    #-------------------------------------------------------------------------------------------
    def plot_jet_raa(self):

        for R in [0.2, 0.4]:
        
            # Get experimental data
            h_data_list = []
            if R==0.2:
                f = ROOT.TFile(self.inclusive_jet_observables['pt_alice']['hepdata_0_10_R02'], 'READ')
                dir = f.Get('Table 30')
            elif R==0.4:
                f = ROOT.TFile(self.inclusive_jet_observables['pt_alice']['hepdata_0_10_R04'], 'READ')
                dir = f.Get('Table 31')
            h_data = dir.Get('Graph1D_y1')
            h_data_list.append([h_data, '0-10%'])
            f.Close()
            
            # Plot
            self.plot_raa(raa_type='jet',
                          hname = f'h_jet_pt_alice_R{R}{self.suffix}Scaled',
                          h_data_list=h_data_list,
                          eta_cut=np.round(self.inclusive_jet_observables['pt_alice']['eta_cut_R']-R, decimals=1),
                          data_centralities=['0-10'],
                          mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                          xtitle="#it{p}_{T,jet} (GeV/#it{c})",
                          ytitle = '#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}} #left[(GeV/c)^{-1}#right]',
                          ymax=1.8,
                          outputfilename=f'h_jet_RAA_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                          R=R,
                          do_chi2=True)
                          
    #-------------------------------------------------------------------------------------------
    def plot_chjet_g(self):
        
        # Get experimental data
        h_data_list = []
        f = ROOT.TFile(self.inclusive_chjet_observables['g_alice']['hepdata'], 'READ')
        dir = f.Get('Table 11')
        h_data = dir.Get('Graph1D_y1')
        h_data_list.append([h_data, '0-10%'])
        f.Close()
        
        # Plot
        R = 0.2
        xtitle="#it{g}"
        self.plot_raa(raa_type='chjet_g',
                      hname = f'h_chjet_g_alice_R{R}{self.suffix}Scaled',
                      h_data_list=None,
                      h_data_list_ratio=h_data_list,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                      xtitle=xtitle,
                      ytitle = f'#frac{{1}}{{#sigma}} #frac{{d#sigma}}{{d#it{{{xtitle}}}}}',
                      ymax=2.8,
                      outputfilename=f'h_chjet_g_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                      R=R,
                      self_normalize=True,
                      do_chi2=True)
                      
    #-------------------------------------------------------------------------------------------
    def plot_chjet_mass(self):
        
        # Get experimental data
        h_data_list = []
        
        f_AA = ROOT.TFile(self.inclusive_chjet_observables['mass_alice']['hepdata_AA'], 'READ')
        dir = f_AA.Get('Table 4')
        h_data_PbPb = dir.Get('Graph1D_y1')
        
        f_pp = ROOT.TFile(self.inclusive_chjet_observables['mass_alice']['hepdata_pp'], 'READ')
        dir = f_pp.Get('Table 1')
        h_data_pp = dir.Get('Graph1D_y1')
                
        h_data_list.append([h_data_PbPb, '0-10%'])
        f_AA.Close()
        f_pp.Close()
        
        # Plot
        R = 0.4
        xtitle="#it{m}"
        self.plot_raa(raa_type='chjet_mass',
                      hname = f'h_chjet_mass_alice_R{R}{self.suffix}Scaled',
                      h_data_list=None,
                      h_data_list_ratio=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                      xtitle=xtitle,
                      ytitle = f'#frac{{1}}{{#sigma}} #frac{{d#sigma}}{{d#it{{{xtitle}}}}}',
                      ymax=2.8,
                      outputfilename=f'h_chjet_mass_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                      R=R,
                      self_normalize=True,
                      do_chi2=True)
                              
    #-------------------------------------------------------------------------------------------
    def plot_chjet_zg(self):
        
        # Plot
        R = 0.4
        xtitle="#it{z}_{g}"
        self.plot_raa(raa_type='chjet_zg',
                      hname = f'h_chjet_zg_alice_R{R}{self.suffix}Scaled',
                      h_data_list=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                      xtitle=xtitle,
                      ytitle = f'#frac{{1}}{{#sigma}} #frac{{d#sigma}}{{d#it{{{xtitle}}}}}',
                      ymax=2.8,
                      outputfilename=f'h_chjet_zg_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                      R=R,
                      self_normalize=True,
                      do_chi2=False)
                      
    #-------------------------------------------------------------------------------------------
    def plot_chjet_tg(self):
        
        # Plot
        R = 0.4
        xtitle="#it{#theta}_{g}"
        self.plot_raa(raa_type='chjet_tg',
                      hname = f'h_chjet_tg_alice_R{R}{self.suffix}Scaled',
                      h_data_list=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                      xtitle=xtitle,
                      ytitle = f'#frac{{1}}{{#sigma}} #frac{{d#sigma}}{{d#it{{{xtitle}}}}}',
                      ymax=2.8,
                      outputfilename=f'h_chjet_tg_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                      R=R,
                      self_normalize=True,
                      do_chi2=False)
                      
    #-------------------------------------------------------------------------------------------
    def plot_semi_inclusive_chjet_IAA(self):
            
        # Get experimental data
        R=0.4
        c_ref = 0.96 # R02: 0.99, R04: 0.96, R05: 0.93
        h_data_list = []
        if R == 0.2:
            f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_IAA_276_R02'], 'READ')
            dir = f.Get('Table 32')
        elif R == 0.4:
            f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_IAA_276_R04'], 'READ')
            dir = f.Get('Table 33')
        elif R == 0.5:
            f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_IAA_276_R05'], 'READ')
            dir = f.Get('Table 34')
        h_data = dir.Get('Graph1D_y1')
        h_data_list.append([h_data, '0-10%'])
        f.Close()
        
        if self.sqrts == 2760:
            hname_ntrigger = f'h_semi_inclusive_chjet_hjet_ntrigger_alice_R{R}{self.suffix}Scaled'
            hname_high = f'h_semi_inclusive_chjet_IAA_highTrigger_alice_R{R}_276{self.suffix}Scaled'
            hname_low = f'h_semi_inclusive_chjet_IAA_lowTrigger_alice_R{R}_276{self.suffix}Scaled'
            
            # Get JETSCAPE pp prediction
            filename_pp = os.path.join(self.output_dir, f'{self.dir_pp}/AnalysisResultsFinal.root')
            f_pp = ROOT.TFile(filename_pp, 'READ')
            h_pp_ntrigger = f_pp.Get(hname_ntrigger)
            h_pp_ntrigger.SetDirectory(0)
            h_pp_high = f_pp.Get(hname_high)
            h_pp_high.SetDirectory(0)
            h_pp_low = f_pp.Get(hname_low)
            h_pp_low.SetDirectory(0)
            f_pp.Close()
            
            # Get JETSCAPE AA prediction
            filename = os.path.join(self.output_dir, f'{self.dir_AA}/AnalysisResultsFinal.root')
            f_AA = ROOT.TFile(filename, 'READ')
            h_AA_ntrigger = f_AA.Get(hname_ntrigger)
            h_AA_ntrigger.SetDirectory(0)
            h_AA_high = f_AA.Get(hname_high)
            h_AA_high.SetDirectory(0)
            h_AA_low = f_AA.Get(hname_low)
            h_AA_low.SetDirectory(0)
            f_AA.Close()
        elif self.sqrts == 5020:
            hname_ntrigger = f'h_semi_inclusive_chjet_hjet_ntrigger_alice_R{R}{self.suffix}Scaled'
            hname_high = f'h_semi_inclusive_chjet_IAA_dphi_highTrigger_alice_R{R}_502{self.suffix}Scaled'
            hname_low = f'h_semi_inclusive_chjet_IAA_dphi_lowTrigger_alice_R{R}_502{self.suffix}Scaled'
            
            # Get JETSCAPE pp prediction
            filename_pp = os.path.join(self.output_dir, f'{self.dir_pp}/AnalysisResultsFinal.root')
            f_pp = ROOT.TFile(filename_pp, 'READ')
            h_pp_ntrigger = f_pp.Get(hname_ntrigger)
            h_pp_ntrigger.SetDirectory(0)
            
            h_pp_high2d = f_pp.Get(hname_high)
            h_pp_high2d.SetDirectory(0)
            h_pp_high2d.GetYaxis().SetRangeUser(np.pi-0.6, np.pi)
            h_pp_high = h_pp_high2d.ProjectionX()
            h_pp_high.Rebin(20)
            h_pp_high.SetDirectory(0)
                        
            h_pp_low2d = f_pp.Get(hname_low)
            h_pp_low2d.SetDirectory(0)
            h_pp_low2d.GetYaxis().SetRangeUser(np.pi-0.6, np.pi)
            h_pp_low = h_pp_low2d.ProjectionX()
            h_pp_low.Rebin(20)
            h_pp_low.SetDirectory(0)
            f_pp.Close()

            # Get JETSCAPE AA prediction
            filename = os.path.join(self.output_dir, f'{self.dir_AA}/AnalysisResultsFinal.root')
            f_AA = ROOT.TFile(filename, 'READ')
            h_AA_ntrigger = f_AA.Get(hname_ntrigger)
            h_AA_ntrigger.SetDirectory(0)
            
            h_AA_high2d = f_AA.Get(hname_high)
            h_AA_high2d.SetDirectory(0)
            h_AA_high2d.GetYaxis().SetRangeUser(np.pi-0.6, np.pi)
            h_AA_high = h_AA_high2d.ProjectionX()
            h_AA_high.Rebin(20)
            h_AA_high.SetDirectory(0)
              
            h_AA_low2d = f_AA.Get(hname_low)
            h_AA_low2d.SetDirectory(0)
            h_AA_low2d.GetYaxis().SetRangeUser(np.pi-0.6, np.pi)
            h_AA_low = h_AA_low2d.ProjectionX()
            h_AA_low.Rebin(20)
            h_AA_low.SetDirectory(0)
            f_AA.Close()
        
        # Delta recoil
        n_trig_high_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(30.))
        n_trig_low_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(8.5))
        n_trig_high_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(30.))
        n_trig_low_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(8.5))
        print(f'n_trig_high_pp: {n_trig_high_pp}')
        print(f'n_trig_low_pp: {n_trig_low_pp}')
        print(f'n_trig_high_AA: {n_trig_high_AA}')
        print(f'n_trig_low_AA: {n_trig_low_AA}')

        h_pp_high.Scale(1./n_trig_high_pp, 'width')
        h_pp_low.Scale(1./n_trig_low_pp, 'width')
        h_AA_high.Scale(1./n_trig_high_AA, 'width')
        h_AA_low.Scale(1./n_trig_low_AA, 'width')
        
        h_delta_recoil_pp = h_pp_high.Clone('h_delta_recoil_pp')
        h_delta_recoil_pp.Add(h_pp_low, -1)
        
        h_delta_recoil_AA = h_AA_high.Clone('h_delta_recoil_AA')
        h_delta_recoil_AA.Add(h_AA_low, -1*c_ref)
 
        # Plot
        xtitle="#it{p}_{T,ch} (GeV/#it{c})"
        self.plot_raa_ratio(raa_type='hjet_IAA',
                            h_pp=h_delta_recoil_pp,
                            h_AA=h_delta_recoil_AA,
                            h_data_list=h_data_list,
                            eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                            data_centralities=['0-10'],
                            mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                            xtitle=xtitle,
                            ytitle = '#Delta_{recoil}',
                            ymax=1.8,
                            outputfilename=f'h_semi_inclusive_chjet_IAA_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                            R=R,
                            do_chi2=True)

    #-------------------------------------------------------------------------------------------
    def plot_semi_inclusive_chjet_dphi(self):
            
        # Get experimental data
        R=0.4
        c_ref = 0.96 # R02: 0.99, R04: 0.96, R05: 0.93
        h_data_list = []
        f = ROOT.TFile(self.semi_inclusive_chjet_observables['hjet_alice']['hepdata_dphi_276'], 'READ')
        dir = f.Get('Table 37')
        h_data = dir.Get('Graph1D_y1')
        h_data_list.append([h_data, 'pp'])
        f.Close()
        
        hname_ntrigger = f'h_semi_inclusive_chjet_hjet_ntrigger_alice_R{R}{self.suffix}Scaled'
        hname_high = f'h_semi_inclusive_chjet_dphi_highTrigger_alice_R{R}_276{self.suffix}Scaled'
        hname_low = f'h_semi_inclusive_chjet_dphi_lowTrigger_alice_R{R}_276{self.suffix}Scaled'
        
        # Get JETSCAPE pp prediction
        filename_pp = os.path.join(self.output_dir, f'{self.dir_pp}/AnalysisResultsFinal.root')
        f_pp = ROOT.TFile(filename_pp, 'READ')
        h_pp_ntrigger = f_pp.Get(hname_ntrigger)
        h_pp_ntrigger.SetDirectory(0)
        h_pp_high = f_pp.Get(hname_high)
        h_pp_high.SetDirectory(0)
        h_pp_low = f_pp.Get(hname_low)
        h_pp_low.SetDirectory(0)
        f_pp.Close()
        
        # Get JETSCAPE AA prediction
        filename = os.path.join(self.output_dir, f'{self.dir_AA}/AnalysisResultsFinal.root')
        f_AA = ROOT.TFile(filename, 'READ')
        h_AA_ntrigger = f_AA.Get(hname_ntrigger)
        h_AA_ntrigger.SetDirectory(0)
        h_AA_high = f_AA.Get(hname_high)
        h_AA_high.SetDirectory(0)
        h_AA_low = f_AA.Get(hname_low)
        h_AA_low.SetDirectory(0)
        f_AA.Close()
        
        # Delta recoil
        n_trig_high_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(30.))
        n_trig_low_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(8.5))
        n_trig_high_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(30.))
        n_trig_low_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(8.5))
        print(f'n_trig_high_pp: {n_trig_high_pp}')
        print(f'n_trig_low_pp: {n_trig_low_pp}')
        print(f'n_trig_high_AA: {n_trig_high_AA}')
        print(f'n_trig_low_AA: {n_trig_low_AA}')

        h_pp_high.Scale(1./n_trig_high_pp, 'width')
        h_pp_low.Scale(1./n_trig_low_pp, 'width')
        h_AA_high.Scale(1./n_trig_high_AA, 'width')
        h_AA_low.Scale(1./n_trig_low_AA, 'width')
        
        h_delta_Phi_pp = h_pp_high.Clone('h_delta_Phi_pp')
        h_delta_Phi_pp.Add(h_pp_low, -1)
        
        h_delta_Phi_AA = h_AA_high.Clone('h_delta_Phi_AA')
        h_delta_Phi_AA.Add(h_AA_low, -1*c_ref)
 
        # Plot
        xtitle="#Delta #it{#varphi}"
        self.plot_raa_ratio(raa_type='hjet_dphi',
                            h_pp=h_delta_Phi_pp,
                            h_AA=h_delta_Phi_AA,
                            h_data_list=h_data_list,
                            eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                            data_centralities=['0-10'],
                            mc_centralities=[f'{self.min_cent}-{self.max_cent}'],
                            xtitle=xtitle,
                            ytitle = '#Phi(#Delta#it{#varphi})',
                            ymax=0.1,
                            outputfilename=f'h_semi_inclusive_chjet_dphi_alice_R{R}_{self.sqrts}_{self.min_cent}-{self.max_cent}{self.file_format}',
                            R=R,
                            do_chi2=True)
     
    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, observable_type, observable, logy = False):
    
        c = ROOT.TCanvas('c', 'c', 600, 650)
        c.Draw()
        c.cd()
        
        # Distribution
        pad2_dy = 0.45
        pad1 = ROOT.TPad('myPad', 'The pad',0,pad2_dy,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(0,1)
        pad1.Draw()
        if logy:
            pad1.SetLogy()
        pad1.cd()
        
        legend = ROOT.TLegend(0.58,0.65,0.75,0.8)
        self.plot_utils.setup_legend(legend, 0.055, sep=-0.1)
            
        self.bins = np.array(self.observable_settings['jetscape_distribution'].GetXaxis().GetXbins())
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, self.bins[0], self.bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.SetMaximum(self.y_max)
        myBlankHisto.SetMinimum(self.y_min) # Don't draw 0 on top panel
        myBlankHisto.GetYaxis().SetTitleSize(0.08)
        myBlankHisto.GetYaxis().SetTitleOffset(1.1)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')

        # Ratio
        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, pad2_dy)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.21)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.04)
        pad2.SetTicks(0,1)
        pad2.Draw()
        pad2.cd()
              
        myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
        myBlankHisto2.SetYTitle('#frac{JETSCAPE}{Data}')
        myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.GetXaxis().SetTitleSize(26)
        myBlankHisto2.GetXaxis().SetTitleFont(43)
        myBlankHisto2.GetXaxis().SetTitleOffset(2.3)
        myBlankHisto2.GetXaxis().SetLabelFont(43)
        myBlankHisto2.GetXaxis().SetLabelSize(22)
        myBlankHisto2.GetYaxis().SetTitleSize(28)
        myBlankHisto2.GetYaxis().SetTitleFont(43)
        myBlankHisto2.GetYaxis().SetTitleOffset(2.)
        myBlankHisto2.GetYaxis().SetLabelFont(43)
        myBlankHisto2.GetYaxis().SetLabelSize(20)
        myBlankHisto2.GetYaxis().SetNdivisions(505)
        myBlankHisto2.GetYaxis().SetRangeUser(self.y_ratio_min, self.y_ratio_max)
        myBlankHisto2.Draw('')
                
        # Draw JETSCAPE
        pad1.cd()
        self.observable_settings['jetscape_distribution'].SetFillColor(self.jetscape_color)
        self.observable_settings['jetscape_distribution'].SetFillColorAlpha(self.jetscape_color, self.alpha)
        self.observable_settings['jetscape_distribution'].SetFillStyle(1001)
        self.observable_settings['jetscape_distribution'].SetMarkerSize(0.)
        self.observable_settings['jetscape_distribution'].SetMarkerStyle(0)
        self.observable_settings['jetscape_distribution'].SetLineWidth(0)
        self.observable_settings['jetscape_distribution'].DrawCopy('E3 same')
        legend.AddEntry(self.observable_settings['jetscape_distribution'], 'JETSCAPE', 'f')

        # Draw distribution
        self.observable_settings['data_distribution'].SetMarkerSize(self.marker_size)
        self.observable_settings['data_distribution'].SetMarkerStyle(self.data_marker)
        self.observable_settings['data_distribution'].SetMarkerColor(self.data_color)
        self.observable_settings['data_distribution'].SetLineStyle(self.line_style)
        self.observable_settings['data_distribution'].SetLineWidth(self.line_width)
        self.observable_settings['data_distribution'].SetLineColor(self.data_color)
        self.observable_settings['data_distribution'].Draw('PE Z X0 same')
        legend.AddEntry(self.observable_settings['data_distribution'], 'Data', 'PE')
        
        legend.Draw()
        
        # Draw ratio
        pad2.cd()
        self.observable_settings['ratio'].SetFillColor(self.jetscape_color)
        self.observable_settings['ratio'].SetFillColorAlpha(self.jetscape_color, self.alpha)
        self.observable_settings['ratio'].SetFillStyle(1001)
        self.observable_settings['ratio'].SetMarkerSize(0.)
        self.observable_settings['ratio'].SetMarkerStyle(0)
        self.observable_settings['ratio'].SetLineWidth(0)
        self.observable_settings['ratio'].Draw('E3 same')

        line = ROOT.TLine(self.bins[0], 1, self.bins[-1], 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        
        pad1.cd()
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        x = 0.35
        text_latex.SetTextSize(0.065)
        #text = f'#bf{{{self.observable}}} #sqrt{{#it{{s_{{#it{{NN}}}}}}}} = {self.sqrts/1000.} TeV'
        text = f'#bf{{{observable_type}_{observable}}} #sqrt{{#it{{s}}}} = {self.sqrts/1000.} TeV'
        text_latex.DrawLatex(x, 0.83, text)

        #text = 'Charged-particle jets'
        #text_latex.DrawLatex(x, 0.75, text)
        
        #text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {0.9-self.jetR}'
        #text_latex.DrawLatex(x, 0.67, text)
        
        #pt_label = None
        #if self.observable in ['girth', 'mass', 'zg', 'theta_g']:
        #    pt_label = f'{self.pt[0]} < #it{{p}}_{{T, ch jet}} < {self.pt[1]} GeV/#it{{c}}'
        #    text_latex.DrawLatex(x, 0.57, pt_label)
            
        #if self.observable in ['zg', 'theta_g']:
        #    grooming_label = f'Soft Drop #it{{z}}_{{cut}}={self.zcut}, #it{{#beta}}=0'
        #    text_latex.DrawLatex(x, 0.47, grooming_label)
        
        c.SaveAs(os.path.join(self.output_dir, f'{self.hname}.pdf'))
        c.Close()

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('Executing plot_results_STAT.py...')
    print('')

    # Define arguments
    parser = argparse.ArgumentParser(description='Plot JETSCAPE events')
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
        default='pp_5020_plot/histograms_5020_merged.root',
        help='Input file'
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

    analysis = PlotResults(config_file=args.configFile, input_file=args.inputFile)
    analysis.plot_results()
