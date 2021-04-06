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

# Base class
sys.path.append('.')
from jetscape_analysis.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotResults(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', output_dir='', **kwargs):
        super(PlotResults, self).__init__(**kwargs)
        self.output_dir = output_dir

        self.data_color = ROOT.kGray+3
        self.data_markers = [21, 20]
        self.markers = [21, 20, 34, 33, 22, 23]
        self.alpha = 0.7
        self.theory_colors = [ROOT.kViolet-8, ROOT.kRed-7, ROOT.kTeal-8,
                              ROOT.kViolet-8, ROOT.kRed-7, ROOT.kTeal-8,
                              ROOT.kBlue-7, ROOT.kBlue-7, ROOT.kBlue-7]
        self.file_format = '.pdf'
        self.debug_level = 0

        # Read config file
        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        #------------------------------------------------------
        # JETSCAPE predictions
        self.jet_R = config['jet_R']
        self.hadron_observables = config['hadron']
        self.inclusive_jet_observables = config['inclusive_jet']
        self.inclusive_chjet_observables = config['inclusive_chjet']
        self.semi_inclusive_chjet_observables = config['semi_inclusive_chjet']
        
        self.dir_pp = '5020_PP'
        self.dir_AA = 'OutputFile_Type5_qhatA10_B100_5020_PbPb_0-10_0.30_2.0_1'
        
        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):

        self.setOptions()
        ROOT.gROOT.ForceStyle()
        
        #----------------------------
        # TG3

        # Charged particle RAA
        self.plot_hadron_raa()

        # Inclusive full jet RAA
        self.plot_jet_raa()
        
        # Charged jet g
        self.plot_chjet_g()
        
        # Charged jet mass
        self.plot_chjet_mass()
        
        # Charged Soft Drop zg
        self.plot_chjet_zg()
        
        # Charged Soft Drop theta_g
        self.plot_chjet_tg()
        
        # h-jet
        self.plot_semi_inclusive_chjet_IAA()
        #self.plot_semi_inclusive_chjet_dphi()
        
        #----------------------------
        # Some extra fun
        have_fun = True
        
        if have_fun:
        
            # Charged jet angularity
            self.plot_chjet_angularity()
            
            # Subjet z
            self.plot_chjet_subjetz()
            
            # Jet axis
            self.plot_chjet_axis()

    #-------------------------------------------------------------------------------------------
    def plot_hadron_raa(self):
    
        # Get experimental data
        h_data_list = []
        f = ROOT.TFile(self.hadron_observables['pt_alice']['hepdata'], 'READ')
        dir = f.Get('Table 8')
        h_data_0_5 = dir.Get('Graph1D_y1')
        h_data_5_10 = dir.Get('Graph1D_y2')
        h_data_list.append([h_data_0_5, '0-5'])
        h_data_list.append([h_data_5_10, '5-10'])
        f.Close()

        # Plot
        self.plot_raa(raa_type='hadron',
                      hname = 'h_hadron_pt_aliceScaled',
                      h_data_list=h_data_list,
                      eta_cut=np.round(self.hadron_observables['pt_alice']['eta_cut'], decimals=1),
                      data_centralities=['0-5', '5-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{p}_{T} (GeV/#it{c})",
                      ymax=1.8,
                      outputfilename=f'h_hadron_RAA_alice{self.file_format}',
                      do_chi2=True)

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
            h_data_list.append([h_data, '0-10'])
            f.Close()
            
            # Plot
            self.plot_raa(raa_type='jet',
                          hname = f'h_jet_pt_alice_R{R}Scaled',
                          h_data_list=h_data_list,
                          eta_cut=np.round(self.inclusive_jet_observables['pt_alice']['eta_cut_R']-R, decimals=1),
                          data_centralities=['0-10'],
                          mc_centralities=['0-10'],
                          xtitle="#it{p}_{T,jet} (GeV/#it{c})",
                          ymax=1.8,
                          outputfilename=f'h_jet_RAA_alice_R{R}{self.file_format}',
                          R=R,
                          do_chi2=True)
                          
    #-------------------------------------------------------------------------------------------
    def plot_chjet_g(self):
        
        # Get experimental data
        h_data_list = []
        f = ROOT.TFile(self.inclusive_chjet_observables['g_alice']['hepdata'], 'READ')
        dir = f.Get('Table 11')
        h_data = dir.Get('Graph1D_y1')
        h_data_list.append([h_data, '0-10'])
        f.Close()
        
        # Plot
        R = 0.2
        self.plot_raa(raa_type='chjet_g',
                      hname = f'h_chjet_g_alice_R{R}Scaled',
                      h_data_list=h_data_list,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{g}",
                      ymax=2.8,
                      outputfilename=f'h_chjet_g_alice_R{R}{self.file_format}',
                      R=R,
                      do_chi2=True)
                      
    #-------------------------------------------------------------------------------------------
    def plot_chjet_mass(self):
        
        # Get experimental data
        h_data_list = []
        f = ROOT.TFile(self.inclusive_chjet_observables['mass_alice']['hepdata'], 'READ')
        dir = f.Get('Table 4')
        h_data = dir.Get('Graph1D_y1')
        h_data_list.append([h_data, '0-10'])
        f.Close()
        
        # Plot
        R = 0.4
        self.plot_raa(raa_type='chjet_mass',
                      hname = f'h_chjet_mass_alice_R{R}Scaled',
                      h_data_list=h_data_list,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{m}",
                      ymax=2.8,
                      outputfilename=f'h_chjet_mass_alice_R{R}{self.file_format}',
                      R=R,
                      do_chi2=True)
                              
    #-------------------------------------------------------------------------------------------
    def plot_chjet_zg(self):
        
        # Plot
        R = 0.2
        self.plot_raa(raa_type='chjet_zg',
                      hname = f'h_chjet_zg_alice_R{R}Scaled',
                      h_data_list=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{z}_{g}",
                      ymax=2.8,
                      outputfilename=f'h_chjet_zg_alice_R{R}{self.file_format}',
                      R=R,
                      do_chi2=False)
                      
    #-------------------------------------------------------------------------------------------
    def plot_chjet_tg(self):
        
        # Plot
        R = 0.2
        self.plot_raa(raa_type='chjet_tg',
                      hname = f'h_chjet_tg_alice_R{R}Scaled',
                      h_data_list=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{#theta}_{g}",
                      ymax=2.8,
                      outputfilename=f'h_chjet_tg_alice_R{R}{self.file_format}',
                      R=R,
                      do_chi2=False)
                      
    #-------------------------------------------------------------------------------------------
    def plot_semi_inclusive_chjet_IAA(self):
    
        hepdata_IAA_276: 'data/TG3/hjet_IAA/HEPData-ins1376027-v1-Table_33.root'
        hepdata_dphi_276: 'data/TG3/hjet_delta_phi/HEPData-ins1376027-v1-Table_37.root'
        
        # Plot
        R = 0.2
        self.plot_raa(raa_type='chjet_tg',
                      hname = f'h_chjet_tg_alice_R{R}Scaled',
                      h_data_list=None,
                      eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                      data_centralities=['0-10'],
                      mc_centralities=['0-10'],
                      xtitle="#it{#theta}_{g}",
                      ymax=2.8,
                      outputfilename=f'h_chjet_tg_alice_R{R}{self.file_format}',
                      R=R,
                      do_chi2=False)
                      
    #-------------------------------------------------------------------------------------------
    def plot_chjet_angularity(self):
        
        for R in [0.2, 0.4]:
            for label in ['groomed', 'ungroomed']:
                for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                
                    self.plot_raa(raa_type='chjet_angularity',
                                  hname = f'h_chjet_angularity_{label}_alice_R{R}_alpha{alpha}Scaled',
                                  h_data_list=None,
                                  eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                                  data_centralities=['0-10'],
                                  mc_centralities=['0-10'],
                                  xtitle=f"#it{{#lambda}}_{{{alpha},{label}}}",
                                  ymax=2.8,
                                  outputfilename=f'h_chjet_angularity_{label}_alice_R{R}_alpha{alpha}{self.file_format}',
                                  R=R,
                                  do_chi2=False)
                                  
    #-------------------------------------------------------------------------------------------
    def plot_chjet_subjetz(self):
        
        for R in [0.2, 0.4]:
            for r in self.inclusive_chjet_observables['subjetz_alice']['r']:
                if r < R:
                
                    self.plot_raa(raa_type='chjet_subjetz',
                                  hname = f'h_chjet_subjetz_alice_R{R}_r{r}Scaled',
                                  h_data_list=None,
                                  eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                                  data_centralities=['0-10'],
                                  mc_centralities=['0-10'],
                                  xtitle=f"#it{{z}}_{{{r}}}",
                                  ymax=2.8,
                                  outputfilename=f'h_chjet_subjetz_alice_R{R}_r{r}{self.file_format}',
                                  R=R,
                                  do_chi2=False)
                                  
    #-------------------------------------------------------------------------------------------
    def plot_chjet_axis(self):
        
        for R in [0.2, 0.4]:

            self.plot_raa(raa_type='chjet_axis',
                          hname = f'h_chjet_axis_Standard_WTA_alice_R{R}Scaled',
                          h_data_list=None,
                          eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                          data_centralities=['0-10'],
                          mc_centralities=['0-10'],
                          xtitle="#Delta#it{R}_{Standard_WTA}",
                          ymax=2.8,
                          outputfilename=f'h_chjet_axis_Standard_WTA_alice_R{R}{self.file_format}',
                          R=R,
                          do_chi2=False)
                          
            self.plot_raa(raa_type='chjet_axis',
                          hname = f'h_chjet_axis_Standard_SD_alice_R{R}Scaled',
                          h_data_list=None,
                          eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                          data_centralities=['0-10'],
                          mc_centralities=['0-10'],
                          xtitle="#Delta#it{R}_{Standard_SD}",
                          ymax=2.8,
                          outputfilename=f'h_chjet_axis_Standard_SD_alice_R{R}{self.file_format}',
                          R=R,
                          do_chi2=False)
                          
            self.plot_raa(raa_type='chjet_axis',
                          hname = f'h_chjet_axis_SD_WTA_alice_R{R}Scaled',
                          h_data_list=None,
                          eta_cut=np.round(self.inclusive_chjet_observables['eta_cut_alice_R']-R, decimals=1),
                          data_centralities=['0-10'],
                          mc_centralities=['0-10'],
                          xtitle="#Delta#it{R}_{SD_WTA}",
                          ymax=2.8,
                          outputfilename=f'h_chjet_axis_SD_WTA_alice_R{R}{self.file_format}',
                          R=R,
                          do_chi2=False)

    #-------------------------------------------------------------------------------------------
    def plot_raa(self, raa_type, hname, h_data_list, eta_cut, data_centralities, mc_centralities,
                 xtitle, ymax, outputfilename,R=None, do_chi2=False):

        # Get JETSCAPE pp prediction
        filename_pp = os.path.join(self.output_dir, f'{self.dir_pp}/AnalysisResultsFinal.root')
        f_pp = ROOT.TFile(filename_pp, 'READ')
        h_pp = f_pp.Get(hname)
        h_pp.SetDirectory(0)
        f_pp.Close()
        
        # Get JETSCAPE AA prediction
        filename = os.path.join(self.output_dir, f'{self.dir_AA}/AnalysisResultsFinal.root')
        f_AA = ROOT.TFile(filename, 'READ')
        h_AA = f_AA.Get(hname)
        h_AA.SetDirectory(0)
        f_AA.Close()

        # For hadrons, impose a 1 GeV minimum, and subtract the recoil hadrons
        if raa_type == 'hadron':
            # Impose 1 Gev minimum
            h_pp_xbins = np.array(h_pp.GetXaxis().GetXbins())
            h_pp_xbins = h_pp_xbins[(h_pp_xbins>=1.)]
            h_pp = h_pp.Rebin(h_pp_xbins.size-1, f'{hname}_pp_rebinned', h_pp_xbins)

            h_AA_xbins = np.array(h_AA.GetXaxis().GetXbins())
            h_AA_xbins = h_AA_xbins[(h_AA_xbins>=1.)]
            h_AA = h_AA.Rebin(h_AA_xbins.size-1, f'{hname}_{self.dir_AA}rebinned', h_AA_xbins)
        
            # Subtract holes
            filename = os.path.join(self.output_dir, f'{self.dir_AA}/AnalysisResultsFinal.root')
            f_AA = ROOT.TFile(filename, 'READ')
            h_recoil_name = 'h_hadron_pt_recoilsScaled'
            h_recoil = f_AA.Get(h_recoil_name)
            h_recoil.SetDirectory(0)
            h_recoil_rebinned = h_recoil.Rebin(h_AA_xbins.size-1, f'{h_recoil_name}_{self.dir_AA}', h_AA_xbins)
            h_AA.Add(h_AA, h_recoil_rebinned, 1, -1)
            h_AA.SetDirectory(0)
            f_AA.Close()

        # Plot the ratio
        if h_pp and h_AA:
            if raa_type == 'hadron':
                output_filename = os.path.join(self.output_dir, f'h_hadron_pt_alice{self.file_format}')
                ytitle = '#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}} #left[(GeV/c)^{-1}#right]'
                h_RAA = self.plot_ratio(h_pp, h_AA, output_filename, xtitle, ytitle, cent=mc_centralities[0],
                                        eta_cut=eta_cut, label=raa_type, logy=True)
            elif raa_type == 'jet':
                output_filename = os.path.join(self.output_dir, f'h_jet_pt_alice_R{R}{self.file_format}')
                ytitle = '#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}} #left[(GeV/c)^{-1}#right]'
                h_RAA = self.plot_ratio(h_pp, h_AA, output_filename, xtitle, ytitle, cent=mc_centralities[0],
                                        eta_cut=eta_cut, label=raa_type, R=R, logy=True)
            elif raa_type in ['chjet_g', 'chjet_mass', 'chjet_zg', 'chjet_tg', 'chjet_angularity', 'chjet_subjetz', 'chjet_axis']:
                output_filename = os.path.join(self.output_dir, f'ratio_{outputfilename}')
                ytitle = f'#frac{{dN}}{{d#it{{{xtitle}}}}}'
                h_RAA = self.plot_ratio(h_pp, h_AA, output_filename, xtitle, ytitle, cent=mc_centralities[0],
                                        eta_cut=eta_cut, label=raa_type, R=R, save_plot = (raa_type in ['chjet_g', 'chjet_mass']))
                if raa_type == 'chjet_mass':
                    return

            h_RAA.SetName(f'{h_RAA.GetName()}_alice')

        else:
            print('h_pp_rebinned or h_AA_rebinned not found!!')
        
        #---------------------------------------
        # Draw RAA
        cname = f'c_{outputfilename}'
        c = ROOT.TCanvas(cname, cname, 600, 450)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(0.12);
        c.SetTopMargin(0.05);
        c.SetBottomMargin(0.15);
        c.cd()

        leg = ROOT.TLegend(0.6,0.75,0.8,0.9)
        self.setupLegend(leg,0.03)

        # Draw experimental data
        if h_data_list:
            for i,h_data_entry in enumerate(h_data_list):
                h_data_entry[0].SetMarkerColor(self.data_color)
                h_data_entry[0].SetLineColor(self.data_color)
                h_data_entry[0].SetMarkerStyle(self.data_markers[i])
                h_data_entry[0].SetMarkerSize(2)
                leg.AddEntry(h_data_entry[0], f'ALICE {h_data_entry[1]}%','Pe')

        # Draw JETSCAPE predictions
        h_RAA.SetNdivisions(505)
        h_RAA.GetXaxis().SetTitleSize(24)
        h_RAA.GetXaxis().SetTitleOffset(1.)
        h_RAA.SetXTitle(xtitle)
        h_RAA.GetYaxis().SetTitleSize(24)
        h_RAA.GetYaxis().SetTitleOffset(1.)
        h_RAA.SetYTitle("#it{R}_{AA}")
        h_RAA.GetYaxis().SetRangeUser(0,ymax)
        h_RAA.Draw('PE same')

        chi2 = None
        if do_chi2:
           chi2 = self.calculate_chi2(h_data_entry[0], h_RAA)

        h_RAA.SetMarkerColorAlpha(self.theory_colors[0], self.alpha)
        h_RAA.SetLineColorAlpha(self.theory_colors[0], self.alpha)
        h_RAA.SetLineWidth(1)
        h_RAA.SetMarkerStyle(self.markers[0])
        if chi2:
            leg.AddEntry(h_RAA,f'{mc_centralities[0]}% (#chi^{{2}}: {chi2})','P')
        else:
            leg.AddEntry(h_RAA,f'{mc_centralities[0]}%','P')

        if h_data_list:
            for h_data_entry in h_data_list:
                h_data_entry[0].Draw('PE same')
        h_RAA.Draw('PE same')

        leg.Draw('same')

        line = ROOT.TLine(h_RAA.GetXaxis().GetXmin(),1,h_RAA.GetXaxis().GetXmax(),1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.9
        dy = 0.05
        x = 0.2
        size = 0.04
        system0 = ROOT.TLatex(x,ymax,'#bf{JETSCAPE}')
        system0.SetNDC()
        system0.SetTextSize(size)
        system0.Draw()

        system1 = ROOT.TLatex(x,ymax-dy,'MATTER+LBT')
        system1.SetNDC()
        system1.SetTextSize(size)
        system1.Draw()

        system2 = ROOT.TLatex(x,ymax-2*dy,'Pb-Pb  #sqrt{#it{s}} = 5.02 TeV')
        system2.SetNDC()
        system2.SetTextSize(size)
        system2.Draw()

        if raa_type == 'hadron':
            system3 = ROOT.TLatex(x,ymax-3*dy, f'Charged particles  |#eta| < {eta_cut}')
        else:
            system3 = ROOT.TLatex(x,ymax-3*dy, f'AKT  #it{{R}} = {R}  |#eta_{{jet}}| < {eta_cut}')
        system3.SetNDC()
        system3.SetTextSize(size)
        system3.Draw()
        
        output_filename = os.path.join(self.output_dir, outputfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Plot ratio h1/h2
    def plot_ratio(self, h_pp, h_AA, outputFilename, xtitle, ytitle, cent, eta_cut, label='jet', logy=False, R=None, save_plot=False):

        # Create canvas
        cname = f'c_{outputFilename}'
        c = ROOT.TCanvas(cname,cname,800,850)
        ROOT.SetOwnership(c, False) # For some reason this is necessary to avoid a segfault...some bug in ROOT or pyroot
                                    # Supposedly fixed in https://github.com/root-project/root/pull/3787
        c.cd()
        pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        if logy:
            pad1.SetLogy()
        pad1.Draw()
        pad1.cd()

        # Set pad and histo arrangement
        myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
        myPad.SetLeftMargin(0.22)
        myPad.SetTopMargin(0.04)
        myPad.SetRightMargin(0.04)
        myPad.SetBottomMargin(0.15)

        # Set spectra styles
        h_AA.SetMarkerSize(3)
        h_AA.SetMarkerStyle(33)
        h_AA.SetMarkerColor(600-6)
        h_AA.SetLineStyle(1)
        h_AA.SetLineWidth(2)
        h_AA.SetLineColor(600-6)

        h_pp.SetMarkerSize(2)
        h_pp.SetMarkerStyle(21)
        h_pp.SetMarkerColor(1)
        h_pp.SetLineStyle(1)
        h_pp.SetLineWidth(2)
        h_pp.SetLineColor(1)

        # Draw spectra
        h_pp.SetXTitle(xtitle)
        h_pp.GetYaxis().SetTitleOffset(2.2)
        h_pp.SetYTitle(ytitle)
        if logy:
            h_pp.SetMaximum(h_pp.GetMaximum()*10.)
            h_pp.SetMinimum(h_pp.GetMinimum()/10.)
        else:
            h_pp.SetMaximum(h_pp.GetMaximum()*2.)
            h_pp.SetMinimum(h_pp.GetMinimum()/2.)

        h_pp.Draw('PE X0 same')
        h_AA.Draw('PE X0 same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # Add legends and text
        # # # # # # # # # # # # # # # # # # # # # # # #
        system = ROOT.TLatex(0.49,0.90,'JETSCAPE')
        system.SetNDC()
        system.SetTextSize(0.044)
        system.Draw()

        system2 = ROOT.TLatex(0.49,0.835,'pp  #sqrt{#it{s}} = 5.02 TeV')
        system2.SetNDC()
        system2.SetTextSize(0.044)
        system2.Draw()

        if 'hadron' in label:
            system3 = ROOT.TLatex(0.49 ,0.765, f'Charged particles | #it{{#eta}} | < {eta_cut}')
            system3.SetNDC()
            system3.SetTextSize(0.044)
            system3.Draw()
        else:
            system3 = ROOT.TLatex(0.49 ,0.765, f'Anti-#it{{k}}_{{T}} #it{{R}} = {R} | #it{{#eta}}_{{jet}}| < {eta_cut}')

        myLegend2pp = ROOT.TLegend(0.45,0.6,0.6,0.68)
        self.setupLegend(myLegend2pp,0.04)
        myLegend2pp.AddEntry(h_AA, f'Pb-Pb {cent}%','P')
        myLegend2pp.AddEntry(h_pp,'pp','Pe')
        myLegend2pp.Draw()

        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        # plot ratio
        hRatio = h_AA.Clone()
        hRatio.SetName('hRatio_{}'.format(cname))
        hRatio.Divide(h_pp)
        hRatio.SetMarkerStyle(21)
        hRatio.SetMarkerSize(2)

        hRatio.GetXaxis().SetTitleSize(30)
        hRatio.GetXaxis().SetTitleFont(43)
        hRatio.GetXaxis().SetTitleOffset(4.)
        hRatio.GetXaxis().SetLabelFont(43)
        hRatio.GetXaxis().SetLabelSize(20)
        hRatio.GetXaxis().SetTitle(xtitle)

        hRatio.GetYaxis().SetTitle('#it{R}_{AA}')
        hRatio.GetYaxis().SetTitleSize(20)
        hRatio.GetYaxis().SetTitleFont(43)
        hRatio.GetYaxis().SetTitleOffset(2.2)
        hRatio.GetYaxis().SetLabelFont(43)
        hRatio.GetYaxis().SetLabelSize(20)
        hRatio.GetYaxis().SetNdivisions(505)

        hRatio.SetMinimum(0.)
        hRatio.SetMaximum(1.49)
        hRatio.Draw('P E')

        line = ROOT.TLine(hRatio.GetXaxis().GetXmin(), 1, hRatio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(4)
        line.Draw()

        if save_plot:
            c.SaveAs(outputFilename)

        return hRatio.Clone('{}_{}'.format(hRatio.GetName(), 'new'))

    # Calculate chi2 for two histograms: in this case h1 is TGraphAsymmErrors, h2 is TH1X
    #-------------------------------------------------------------------------------------------
    def calculate_chi2(self, h1, h2):
       count = 0
       chi2 = 0

       # Loop over bins in h1, check if corresponding point exist in h2
       for i in range(h1.GetN()):
          # x, y, exl, exh, eyl, eyh = ROOT.Double(-1), ROOT.Double(-1), -1, -1, -1, -1
          x, y, exl, exh, eyl, eyh = ctypes.c_double(-1), ctypes.c_double(-1), -1, -1, -1, -1
          h1.GetPoint(i, x, y)
          exh = h1.GetErrorXhigh(i)
          exl = h1.GetErrorXlow(i)
          eyh = h1.GetErrorYhigh(i)
          eyl = h1.GetErrorYlow(i)

          # Use average error for now
          ey = (eyl + eyh) / 2

          # Corresponding bin in h2
          b = h2.FindBin(x.value)
          if b <= 0 or b > h2.GetNbinsX():
             continue

          y2 = h2.GetBinContent(b)
          ey2 = h2.GetBinError(b)

          chi2 = chi2 + (y.value - y2)**2 / (ey * ey + ey2 * ey2)
          count = count + 1

       return '{:3.1f}/{:d}'.format(chi2, count)


    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def get_data_from_txt(self, filename, output_filename):

        g = ROOT.TGraphAsymmErrors(filename, "%lg %lg %lg %lg %lg %lg")

        if self.debug_level > 0:
            cname = 'c_{}'.format(filename)
            c = ROOT.TCanvas(cname,cname,800,850)
            g.Draw()
            c.SaveAs(output_filename)
        return g

    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def setupLegend(self, leg, textSize):

        leg.SetTextFont(42);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetFillColor(0);
        leg.SetMargin(0.25);
        leg.SetTextSize(textSize);
        leg.SetEntrySeparation(0.5);

    #-------------------------------------------------------------------------------------------
    def setOptions(self):

        font = 42

        ROOT.gStyle.SetFrameBorderMode(0)
        ROOT.gStyle.SetFrameFillColor(0)
        ROOT.gStyle.SetCanvasBorderMode(0)
        ROOT.gStyle.SetPadBorderMode(0)
        ROOT.gStyle.SetPadColor(10)
        ROOT.gStyle.SetCanvasColor(10)
        ROOT.gStyle.SetTitleFillColor(10)
        ROOT.gStyle.SetTitleBorderSize(1)
        ROOT.gStyle.SetStatColor(10)
        ROOT.gStyle.SetStatBorderSize(1)
        ROOT.gStyle.SetLegendBorderSize(1)

        ROOT.gStyle.SetDrawBorder(0)
        ROOT.gStyle.SetTextFont(font)
        ROOT.gStyle.SetStatFont(font)
        ROOT.gStyle.SetStatFontSize(0.05)
        ROOT.gStyle.SetStatX(0.97)
        ROOT.gStyle.SetStatY(0.98)
        ROOT.gStyle.SetStatH(0.03)
        ROOT.gStyle.SetStatW(0.3)
        ROOT.gStyle.SetTickLength(0.02,"y")
        ROOT.gStyle.SetEndErrorSize(3)
        ROOT.gStyle.SetLabelSize(0.05,"xyz")
        ROOT.gStyle.SetLabelFont(font,"xyz")
        ROOT.gStyle.SetLabelOffset(0.01,"xyz")
        ROOT.gStyle.SetTitleFont(font,"xyz")
        ROOT.gStyle.SetTitleOffset(1.2,"xyz")
        ROOT.gStyle.SetTitleSize(0.045,"xyz")
        ROOT.gStyle.SetMarkerSize(1)
        ROOT.gStyle.SetPalette(1)

        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('Executing plot_results_TG3.py...')
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

    analysis = PlotResults(config_file=args.configFile, output_dir=args.outputDir)
    analysis.plot_results()