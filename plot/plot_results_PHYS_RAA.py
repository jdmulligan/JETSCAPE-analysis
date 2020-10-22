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
        
        # Read config file
        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)
            
        #------------------------------------------------------
        # JETSCAPE predictions
        self.charged_particle_eta_cut = config['charged_particle_eta_cut']
        self.min_hadron_pt = config['min_hadron_pt']
        self.subtract_hadron_recoil = config['subtract_hadron_recoil']
        self.jetR_list = config['jetR']
        self.jet_eta_cut_04 = config['jet_eta_cut_04']
        self.jet_eta_cut_02 = config['jet_eta_cut_02']
        self.debug_level = config['debug_level']
        self.file_format = '.pdf'
        
        self.predictions = {'pp': [], 'central': [], 'semicentral': [],}
        self.subdirs = [x for _,x,_ in os.walk(self.output_dir)][0]
        for subdir in self.subdirs:
        
            if 'PP' in subdir:
                self.predictions['pp'].append(subdir)
                continue
            else:
                tags = str(subdir).split('_')
                centrality = tags[2]
                alpha_s = float(tags[3].replace('R','.'))
                Q_switch = float(tags[4].replace('R','.'))
                
                prediction = [subdir, centrality, alpha_s, Q_switch]
            
                if '0-10' in subdir:
                    self.predictions['central'].append(prediction)
                elif '30-40' in subdir or '40-50' in subdir:
                    self.predictions['semicentral'].append(prediction)
           
        #------------------------------------------------------
        # Charged particle data
        self.file_CMS_hadron_0_5 = config['CMS_hadron_0_5']
        self.file_CMS_hadron_5_10 = config['CMS_hadron_5_10']
        self.file_CMS_hadron_30_50 = config['CMS_hadron_30_50']
        self.file_ATLAS_hadron_0_5 = config['ATLAS_hadron_0_5']
        self.file_ATLAS_hadron_30_40 = config['ATLAS_hadron_30_40']
        self.file_ALICE_hadron = config['ALICE_hadron']

        #------------------------------------------------------
        # Jet data
        self.file_CMS_jet_0_10_R02 = config['CMS_jet_0_10_R02']
        self.file_CMS_jet_0_10_R04 = config['CMS_jet_0_10_R04']
        self.file_CMS_jet_30_50_R02 = config['CMS_jet_30_50_R02']
        self.file_CMS_jet_30_50_R04 = config['CMS_jet_30_50_R04']
        self.file_ATLAS_jet_0_10 = config['ATLAS_jet_0_10']
        self.file_ATLAS_jet_30_40 = config['ATLAS_jet_30_40']
        self.file_ATLAS_jet_40_50 = config['ATLAS_jet_40_50']
        self.file_ALICE_jet_0_10_R02 = config['ALICE_jet_0_10_R02']
        self.file_ALICE_jet_0_10_R04 = config['ALICE_jet_0_10_R04']
        self.file_ALICE_jet_0_10_R02_biasratio = config['ALICE_jet_0_10_R02_biasratio']
        
        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):
      
        self.setOptions()
        ROOT.gROOT.ForceStyle()
        
        # Charged particle histograms
        plot_hadron_histograms = True
        if plot_hadron_histograms:
            self.plot_hadron_histograms()
            
        # Jet histograms
        plot_jet_histograms = True
        if plot_jet_histograms:
            self.plot_jet_histograms()
            
        # Plot some additional histograms
        plot_qa_histograms = True
        if plot_qa_histograms:
            self.plot_qa_histograms()
            
            # Plot ALICE bias ratio in pp
            self.plot_pt_lead_ratio(R=0.2)
            self.plot_pt_lead_ratio(R=0.4)
            
    #-------------------------------------------------------------------------------------------
    def plot_hadron_histograms(self):
    
        # Create multi-panel canvas
        cname = 'c_hadron'
        c = ROOT.TCanvas(cname,cname,1700,900)
        c.cd()
        c.Divide(3, 2)
        
        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []
    
        self.plot_raa('ALICE', c, pad=1, raa_type='hadron', cent_type='central',
                      eta_cut=self.charged_particle_eta_cut[2],
                      data_centralities=['0-5', '5-10'],
                      mc_centralities=['0-10'])
        self.plot_raa('ALICE', c, pad=4, raa_type='hadron', cent_type='semicentral',
                      eta_cut=self.charged_particle_eta_cut[2],
                      data_centralities=['30-40'],
                      mc_centralities=['30-40'])
        self.plot_raa('ATLAS', c, pad=2, raa_type='hadron', cent_type='central',
                      eta_cut=self.charged_particle_eta_cut[1],
                      data_centralities=['0-5'],
                      mc_centralities=['0-10'])
        self.plot_raa('ATLAS', c, pad=5, raa_type='hadron',  cent_type='semicentral',
                      eta_cut=self.charged_particle_eta_cut[1],
                      data_centralities=['30-40'],
                      mc_centralities=['30-40'])
        self.plot_raa('CMS', c, pad=3, raa_type='hadron', cent_type='central',
                      eta_cut=self.charged_particle_eta_cut[0],
                      data_centralities=['0-5', '5-10'],
                      mc_centralities=['0-10'])
        self.plot_raa('CMS', c, pad=6, raa_type='hadron', cent_type='semicentral',
                      eta_cut=self.charged_particle_eta_cut[0],
                      data_centralities=['30-50'],
                      mc_centralities=['30-40', '40-50'])
                    
        output_filename = os.path.join(self.output_dir, 'hHadronRAA{}'.format(self.file_format))
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    def plot_jet_histograms(self):
    
        # Create multi-panel canvas
        cname = 'c_jet'
        c = ROOT.TCanvas(cname,cname,2300,900)
        c.cd()
        c.Divide(4, 2)
        
        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []

        self.plot_raa('ALICE', c, pad=1, raa_type='jet', cent_type='central',
                      eta_cut=self.jet_eta_cut_02[2],
                      data_centralities=['0-10'], mc_centralities=['0-10'], R=0.2)
        self.plot_raa('ALICE', c, pad=5, raa_type='jet', cent_type='central',
                      eta_cut=self.jet_eta_cut_04[2],
                      data_centralities=['0-10'], mc_centralities=['0-10'], R=0.4)
        self.plot_raa('ATLAS', c, pad=2, raa_type='jet', cent_type='central',
                      eta_cut=self.jet_eta_cut_04[1],
                      data_centralities=['0-10'], mc_centralities=['0-10'], R=0.4)
        self.plot_raa('ATLAS', c, pad=6, raa_type='jet', cent_type='semicentral',
                      eta_cut=self.jet_eta_cut_04[1],
                      data_centralities=['30-40'], mc_centralities=['30-40'], R=0.4)
        self.plot_raa('CMS', c, pad=3, raa_type='jet', cent_type='central',
                      eta_cut=self.jet_eta_cut_02[0],
                      data_centralities=['0-10'], mc_centralities=['0-10'], R=0.2)
        self.plot_raa('CMS', c, pad=7, raa_type='jet', cent_type='semicentral',
                      eta_cut=self.jet_eta_cut_02[0],
                      data_centralities=['30-50'], mc_centralities=['30-40'], R=0.2)
        self.plot_raa('CMS', c, pad=4, raa_type='jet', cent_type='central',
                      eta_cut=self.jet_eta_cut_04[0],
                      data_centralities=['0-10'], mc_centralities=['0-10'], R=0.4)
        self.plot_raa('CMS', c, pad=8, raa_type='jet', cent_type='semicentral',
                      eta_cut=self.jet_eta_cut_04[0],
                      data_centralities=['30-50'], mc_centralities=['30-40'], R=0.4)
        
        output_filename = os.path.join(self.output_dir, 'hJetRAA{}'.format(self.file_format))
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    def plot_raa(self, experiment, c, pad, raa_type, cent_type,
                 eta_cut, data_centralities, mc_centralities, R=None):

        # Get JETSCAPE prediction
        if raa_type == 'hadron':
            hname = 'hChargedPt_{}Scaled'.format(experiment)
        elif raa_type == 'jet':
            if experiment in ['ALICE', 'CMS']:
                hname = 'hJetPt_{}_R{}Scaled'.format(experiment, R)
            elif experiment == 'ATLAS':
                if cent_type == 'central':
                    hname = 'hJetPt_{}_binning0_R{}Scaled'.format(experiment, R)
                elif cent_type == 'semicentral':
                    hname = 'hJetPt_{}_binning1_R{}Scaled'.format(experiment, R)

        # pp
        filename_pp = os.path.join(self.output_dir, '{}/AnalysisResultsFinal.root'.format(self.predictions['pp'][0]))
        f_pp = ROOT.TFile(filename_pp, 'READ')
        h_pp = f_pp.Get(hname)
        h_pp.SetDirectory(0)
        f_pp.Close()

        # Impose 1 Gev minimum
        h_pp_xbins = np.array(h_pp.GetXaxis().GetXbins())
        h_pp_xbins = h_pp_xbins[(h_pp_xbins>=self.min_hadron_pt)]
        h_pp_rebinned = h_pp.Rebin(h_pp_xbins.size-1, '{}_pp_rebinned'.format(hname), h_pp_xbins)
        
        # AA
        predictions = sorted(self.predictions[cent_type], key=operator.itemgetter(2))
        if predictions[0][2] == 0.2:
            predictions.append(predictions.pop(0))
            predictions.append(predictions.pop(-2))
        if predictions[1][3] == 3.0:
            predictions_temp1 = predictions[1]
            predictions_temp2 = predictions[2]
            predictions[1] = predictions_temp2
            predictions[2] = predictions_temp1
        
        predictions_to_plot = []
        h_AA = None
        self.h_RAA_list = []
        
        for i,prediction in enumerate(predictions):
            mc_cent = prediction[1]
            alpha_s = prediction[2]
            Q_switch = prediction[3]
            
            if mc_cent in mc_centralities:
                predictions_to_plot.append(prediction)
            else:
                continue
            
            filename = os.path.join(self.output_dir, '{}/AnalysisResultsFinal.root'.format(prediction[0]))
            f_AA = ROOT.TFile(filename, 'READ')
            h_AA = f_AA.Get(hname)
            h_AA.SetDirectory(0)
            
            # Impose 1 Gev minimum
            h_AA_xbins = np.array(h_AA.GetXaxis().GetXbins())
            h_AA_xbins = h_AA_xbins[(h_AA_xbins>=self.min_hadron_pt)]
            h_AA_rebinned = h_AA.Rebin(h_AA_xbins.size-1, '{}_{}rebinned'.format(hname, prediction[0]), h_AA_xbins)
            
            # For hadrons, subtract the recoil hadrons
            if raa_type == 'hadron':
                h_recoil_name = 'hChargedPt_RecoilsScaled'
                h_recoil = f_AA.Get(h_recoil_name)
                h_recoil.SetDirectory(0)
                h_recoil_rebinned = h_recoil.Rebin(h_AA_xbins.size-1, '{}_{}'.format(h_recoil_name, prediction[0]), h_AA_xbins)
                if self.subtract_hadron_recoil:
                    h_AA_rebinned.Add(h_AA_rebinned, h_recoil_rebinned, 1, -1)
            h_AA_rebinned.SetDirectory(0)
            f_AA.Close()

            # Plot the ratio
            if h_pp_rebinned and h_AA_rebinned:
                if raa_type == 'hadron':
                    output_filename = os.path.join(self.output_dir, 'hChHadron_{}_{}_{}{}'.format(cent_type, experiment, i, self.file_format))
                    xtitle = '#it{p}_{T} (GeV/#it{c})'
                    ytitle = '#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}} #left[(GeV/c)^{-1}#right]'
                    h_RAA = self.plot_ratio(h_pp_rebinned, h_AA_rebinned, output_filename, xtitle, ytitle, eta_cut=eta_cut,
                                            cent=mc_cent, alpha_s=alpha_s, Q_switch=Q_switch, label='Hadron')
                elif raa_type == 'jet':
                    output_filename = os.path.join(self.output_dir, 'hJetPt_{}_{}_{}{}'.format(cent_type, experiment, i, self.file_format))
                    xtitle = '#it{p}_{T} (GeV/#it{c})'
                    ytitle = '#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}} #left[(GeV/c)^{-1}#right]'
                    h_RAA = self.plot_ratio(h_pp_rebinned, h_AA_rebinned, output_filename, xtitle, ytitle, eta_cut=eta_cut,
                                            cent=mc_cent, alpha_s=alpha_s, Q_switch=Q_switch, label='Jet', R=R)

                h_RAA.SetName('{}_{}_{}'.format(h_RAA.GetName(), experiment, pad))
                self.plot_list.append(h_RAA)
                self.h_RAA_list.append(h_RAA)
            else:
                print('h_pp_rebinned or h_AA_rebinned not found!!')
                
        if raa_type == 'hadron':
            h_data_list = self.get_hadron_data(experiment, mc_cent)
        elif raa_type == 'jet':
            h_data_list = self.get_jet_data(experiment, mc_cent, R)
        # Plot RAA overlay
        if len(self.h_RAA_list) > 0:
            self.plot_RAA_overlay(raa_type, c, pad, predictions_to_plot, h_data_list, experiment, eta_cut, R)
            
        self.plot_list.append(predictions_to_plot)
        self.plot_list.append(h_data_list)

    #-------------------------------------------------------------------------------------------
    def plot_RAA_overlay(self, raa_type, c, pad, predictions, h_data_list, experiment, eta_cut, R):

        # Create canvas
        c.cd(pad)
        
        # Set pad and histo arrangement
        myPad = ROOT.TPad("myPad_{}".format(raa_type), "The pad{}".format(raa_type),0,0,1,1)
        self.plot_list.append(myPad)

        if raa_type == 'hadron':
            myPad.SetLeftMargin(0.12)
            myPad.SetRightMargin(0.03)
            myPad.SetTopMargin(0.01)
            myPad.SetBottomMargin(0.15)
        elif raa_type == 'jet':
            myPad.SetLeftMargin(0.12)
            myPad.SetRightMargin(0.03)
            myPad.SetTopMargin(0.01)
            myPad.SetBottomMargin(0.15)

        myPad.SetTicks(0,1)
        myPad.Draw()
        myPad.cd()
        
        n_predictions = len(predictions) + len(h_data_list)
        if raa_type == 'hadron':
            leg = ROOT.TLegend(0.55,0.93-0.042*n_predictions,0.75,0.96)
            self.setupLegend(leg,0.04)
        elif raa_type == 'jet':
            leg = ROOT.TLegend(0.53,0.93-0.042*n_predictions,0.73,0.96)
            self.setupLegend(leg,0.04)
        self.plot_list.append(leg)
        
        # Draw experimental data
        for i,h_data_entry in enumerate(h_data_list):
            h_data_entry[0].SetMarkerColor(self.data_color)
            h_data_entry[0].SetLineColor(self.data_color)
            h_data_entry[0].SetMarkerStyle(self.data_markers[i])
            h_data_entry[0].SetMarkerSize(2)
            leg.AddEntry(h_data_entry[0],'{} {}%'.format(experiment, h_data_entry[1]),'Pe')

        # Draw JETSCAPE predictions
        for i,prediction in enumerate(predictions):
            cent = prediction[1]
            alpha_s = prediction[2]
            Q_switch = prediction[3]
            
            if i == 0:
            
                self.h_RAA_list[i].SetNdivisions(505)
                self.h_RAA_list[i].GetXaxis().SetTitleSize(24)
                self.h_RAA_list[i].GetXaxis().SetTitleOffset(2.6)
                self.h_RAA_list[i].SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
                self.h_RAA_list[i].GetYaxis().SetTitleSize(24)
                self.h_RAA_list[i].GetYaxis().SetTitleOffset(1.8)
                self.h_RAA_list[i].SetYTitle("#it{R}_{AA}")
                self.h_RAA_list[i].GetYaxis().SetRangeUser(0,1.8)
                self.h_RAA_list[i].Draw('PE same')

            #self.h_RAA_list[i].SetMarkerColor(self.theory_colors[i])
            self.h_RAA_list[i].SetMarkerColorAlpha(self.theory_colors[i], self.alpha)
            #self.h_RAA_list[i].SetLineColor(self.theory_colors[i])
            self.h_RAA_list[i].SetLineColorAlpha(self.theory_colors[i], self.alpha)
            self.h_RAA_list[i].SetLineWidth(1)
            self.h_RAA_list[i].SetMarkerStyle(self.markers[i])
            leg.AddEntry(self.h_RAA_list[i],'#alpha_{{s}} = {}, #it{{Q}}_{{0}} = {}, {}%'.format(alpha_s, Q_switch, cent),'P')

        for h_data_entry in h_data_list:
            h_data_entry[0].Draw('PE same')
        for i,prediction in enumerate(predictions):
            self.h_RAA_list[i].Draw('PE same')
        
        leg.Draw('same')
        
        line = ROOT.TLine(self.h_RAA_list[i].GetXaxis().GetXmin(),1,self.h_RAA_list[i].GetXaxis().GetXmax(),1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)
        
        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.93
        dy = 0.05
        x = 0.16
        system0 = ROOT.TLatex(x,ymax,'#bf{JETSCAPE}')
        system0.SetNDC()
        system0.SetTextSize(0.044)
        system0.Draw()
        
        system1 = ROOT.TLatex(x,ymax-dy,'MATTER+LBT')
        system1.SetNDC()
        system1.SetTextSize(0.044)
        system1.Draw()
        
        system2 = ROOT.TLatex(x,ymax-2*dy,'Pb-Pb  #sqrt{#it{s}} = 5.02 TeV')
        system2.SetNDC()
        system2.SetTextSize(0.044)
        system2.Draw()
        
        if raa_type == 'hadron':
            system3 = ROOT.TLatex(x,ymax-3*dy, 'Charged particles  |#eta| < {}'.format(eta_cut))
        elif raa_type == 'jet':
            if experiment == 'ATLAS':
                system3 = ROOT.TLatex(x,ymax-3*dy, 'AKT  #it{{R}} = {}  |#it{{y}}_{{jet}}| < {}'.format(R, eta_cut))
            else:
                system3 = ROOT.TLatex(x,ymax-3*dy, 'AKT  #it{{R}} = {}  |#eta_{{jet}}| < {}'.format(R, eta_cut))
        system3.SetNDC()
        system3.SetTextSize(0.044)
        system3.Draw()
        
        self.plot_list.append(system0)
        self.plot_list.append(system1)
        self.plot_list.append(system2)
        self.plot_list.append(system3)

    #-------------------------------------------------------------------------------------------
    def plot_qa_histograms(self):

        # Plot unscaled number of events per pt-hat bin
        self.plot_hist_loop(hname = 'hNevents')
        
        # Plot hadron recoil pt distribution
        self.plot_hist_loop(hname = 'hChargedPt_RecoilsScaled')
        
    #-------------------------------------------------------------------------------------------
    def plot_hist_loop(self, hname):

        cname = 'c_{}'.format(hname)
        c = ROOT.TCanvas(cname,cname,600, 450)
        c.cd()
            
        # Set pad and histo arrangement
        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
        myPad.SetTicks(0,1)
        myPad.SetLeftMargin(0.15)
        myPad.SetRightMargin(0.05)
        myPad.SetBottomMargin(0.12)
        myPad.SetTopMargin(0.05)
        myPad.Draw()
        if 'hChargedPt_RecoilsScaled' in hname:
            myPad.SetLogy()
        myPad.cd()
        
        if 'hChargedPt_RecoilsScaled' in hname:
            leg = ROOT.TLegend(0.4,0.56,0.8,0.92)
            self.setupLegend(leg,0.04)
        else:
            leg = ROOT.TLegend(0.25,0.56,0.65,0.92)
            self.setupLegend(leg,0.04)

        i=0
        h_list = []
        for key,cent_group in self.predictions.items():

            for prediction in cent_group:
                            
                if key == 'pp':
                    dir = prediction
                else:
                    dir = prediction[0]
                
                filename = os.path.join(self.output_dir, '{}/AnalysisResultsFinal.root'.format(dir))
                f = ROOT.TFile(filename, 'READ')
                h = f.Get(hname)
                h.SetDirectory(0)
                h.SetName('{}_{}'.format(h.GetName(), dir))
                h_list.append(h)
                f.Close()
                
                if 'Nevents' in hname:
                    h.GetYaxis().SetTitle('N events (unscaled)')
                    h.GetXaxis().SetTitle('pt hat index')
                    h.GetXaxis().SetTitleOffset(1.)
                    h.GetYaxis().SetRangeUser(-1, 250000)
                    n_event_avg = int(h.Integral()/66.)
                    print('Avg n_events: {} ({})'.format(n_event_avg, dir))
                    if n_event_avg == 120000:
                        h.SetMarkerStyle(21)
                    else:
                        h.SetMarkerStyle(20)
                    leg.AddEntry(h, '{} (avg: {})'.format(dir, n_event_avg))
            
                elif 'hChargedPt_RecoilsScaled' in hname:
                    h.GetYaxis().SetTitle('#frac{dN}{d#it{p}_{T}} (GeV/#it{c})^{-1}')
                    h.GetYaxis().SetTitleOffset(2.8)
                    h.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
                    h.GetYaxis().SetTitleOffset(1.2)
                    h.GetXaxis().SetRangeUser(0,50)
                    h.GetYaxis().SetRangeUser(1e-5,1e12)
                    h.SetMarkerStyle(21)
                    leg.AddEntry(h, dir)

                    #rebin_value = 5
                    #h.Rebin(rebin_value)
                    #h.Scale(1./rebin_value)

                h.SetMarkerColor(ROOT.kBlue+2-i)
                h.Draw('P same')

                i += 1
                
        leg.Draw('same')
        
        output_filename = os.path.join(self.output_dir, '{}{}'.format(hname, self.file_format))
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    def plot_pt_lead_ratio(self, R):

        cname = 'c_ptlead_{}'.format(R)
        c = ROOT.TCanvas(cname,cname,600, 450)
        c.cd()
        
        # Set pad and histo arrangement
        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
        myPad.SetTicks(0,1)
        myPad.SetLeftMargin(0.2)
        myPad.SetRightMargin(0.05)
        myPad.SetBottomMargin(0.12)
        myPad.SetTopMargin(0.05)
        myPad.Draw()
        myPad.cd()
        
        leg = ROOT.TLegend(0.3,0.56,0.7,0.92)
        self.setupLegend(leg,0.04)

        i=0
        h_list = []
        for key,cent_group in self.predictions.items():

            for prediction in cent_group:
                            
                if key == 'pp':
                    dir = prediction
                else:
                    dir = prediction[0]
                
                filename = os.path.join(self.output_dir, '{}/AnalysisResultsFinal.root'.format(dir))
                f = ROOT.TFile(filename, 'READ')
                
                h1name = 'hJetPt_ALICE_R{}Scaled'.format(R)
                h1 = f.Get(h1name)
                h1.SetDirectory(0)
                h1.SetName('{}_{}'.format(h1.GetName(), dir))
  
                h2name = 'hJetPt_ALICE_no_ptlead_cut_R{}Scaled'.format(R)
                h2 = f.Get(h2name)
                h2.SetDirectory(0)
                h2.SetName('{}_{}'.format(h2.GetName(), dir))
  
                f.Close()
 
                hratio = h1.Clone('{}_ratio'.format(h1.GetName()))
                hratio.Divide(h2)
                h_list.append(hratio)
                hratio.GetYaxis().SetRangeUser(0., 2.)
                if R == 0.2:
                    hratio.GetYaxis().SetTitle('bias ratio:   #frac{#it{p}_{T,ch lead} > 5 GeV/#it{c}}{#it{p}_{T,ch lead} > 0 GeV/#it{c}}')
                else:
                    h1.GetYaxis().SetTitle('bias ratio')
                hratio.GetYaxis().SetTitleOffset(2.)
                hratio.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
                hratio.GetXaxis().SetTitleOffset(1.2)
                if key == 'pp':
                    hratio.SetMarkerStyle(20)
                else:
                    hratio.SetMarkerStyle(33)
                hratio.SetMarkerColor(ROOT.kBlue+2-i)
                hratio.Draw('P same')
   
                leg.AddEntry(hratio, dir)

                i += 1
        
        # Draw exp data, for R=0.2 case
        if R == 0.2:
            f = ROOT.TFile(self.file_ALICE_jet_0_10_R02_biasratio, 'READ')
            dir = f.Get('Table 26')
            h_data = dir.Get('Graph1D_y1')
            h_data.SetMarkerStyle(21)
            h_data.SetMarkerColor(self.data_color)
            h_data.SetLineColor(self.data_color)
            h_data.Draw('EP same')
            leg.AddEntry(h_data, 'ALICE pp R=0.2')

        leg.Draw('same')
        
        line = ROOT.TLine(h1.GetXaxis().GetXmin(),1,h1.GetXaxis().GetXmax(),1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        
        output_filename = os.path.join(self.output_dir, 'hBiasRatio_R{}{}'.format(R, self.file_format))
        c.SaveAs(output_filename)
         
    #-------------------------------------------------------------------------------------------
    def get_hadron_data(self, experiment, mc_cent):

        # Get experimental data
        h_data_list = []
        if experiment == 'ALICE':
            f = ROOT.TFile(self.file_ALICE_hadron, 'READ')
            dir = f.Get('Table 8')
            if mc_cent == '0-10':
                h_data_0_5 = dir.Get('Graph1D_y1')
                h_data_5_10 = dir.Get('Graph1D_y2')
                h_data_list.append([h_data_0_5, '0-5'])
                h_data_list.append([h_data_5_10, '5-10'])
            elif mc_cent == '30-40':
                h_data = dir.Get('Graph1D_y5')
                h_data_list.append([h_data, mc_cent])
            f.Close()
        elif experiment == 'ATLAS':
            if mc_cent == '0-10':
                h_data= self.get_data_from_txt(self.file_ATLAS_hadron_0_5,
                                               os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
            elif mc_cent == '30-40':
                h_data= self.get_data_from_txt(self.file_ATLAS_hadron_30_40,
                                               os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
            h_data_list.append([h_data, mc_cent])
        elif experiment == 'CMS':
            if mc_cent == '0-10':
                f_0_5 = ROOT.TFile(self.file_CMS_hadron_0_5, 'READ')
                dir_0_5 = f_0_5.Get('Table 8')
                h_data_0_5 = dir_0_5.Get('Graph1D_y1')
                f_5_10 = ROOT.TFile(self.file_CMS_hadron_5_10, 'READ')
                dir_5_10 = f_5_10.Get('Table 9')
                h_data_5_10 = dir_5_10.Get('Graph1D_y1')
                h_data_list.append([h_data_0_5, '0-5'])
                h_data_list.append([h_data_5_10, '5-10'])
            elif mc_cent in ['30-40', '40-50']:
                f = ROOT.TFile(self.file_CMS_hadron_30_50, 'READ')
                dir = f.Get('Table 11')
                h_data = dir.Get('Graph1D_y1')
                h_data_list.append([h_data, '30-50'])
                f.Close()

        return h_data_list

    #-------------------------------------------------------------------------------------------
    def get_jet_data(self, experiment, mc_cent, R):
    
        # Get experimental data
        h_data_list = []
        if experiment == 'ALICE':
            if mc_cent == '0-10':
                if R==0.2:
                    f = ROOT.TFile(self.file_ALICE_jet_0_10_R02, 'READ')
                    dir = f.Get('Table 30')
                elif R==0.4:
                    f = ROOT.TFile(self.file_ALICE_jet_0_10_R04, 'READ')
                    dir = f.Get('Table 31')
                h_data = dir.Get('Graph1D_y1')
                h_data_list.append([h_data, '0-10'])
            f.Close()
        elif experiment == 'ATLAS':
            if mc_cent == '0-10':
                f = ROOT.TFile(self.file_ATLAS_jet_0_10, 'READ')
                dir = f.Get('Table 19')
                h_data = dir.Get('Graph1D_y1')
                h_data_list.append([h_data, '0-10'])
            elif mc_cent == '30-40':
                f = ROOT.TFile(self.file_ATLAS_jet_30_40, 'READ')
                dir = f.Get('Table 22')
                h_data = dir.Get('Graph1D_y1')
                h_data_list.append([h_data, '30-40'])
            f.Close()
        elif experiment == 'CMS':
            if mc_cent == '0-10':
                if R==0.2:
                    h_data = self.get_data_from_txt(self.file_CMS_jet_0_10_R02,
                                                    os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
                elif R==0.4:
                    h_data = self.get_data_from_txt(self.file_CMS_jet_0_10_R04,
                                                    os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
            elif mc_cent in ['30-40', '40-50']:
                if R==0.2:
                    h_data = self.get_data_from_txt(self.file_CMS_jet_30_50_R02,
                                                    os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
                elif R==0.4:
                    h_data = self.get_data_from_txt(self.file_CMS_jet_30_50_R04,
                                                    os.path.join(self.output_dir, '{}_{}{}'.format(experiment, mc_cent, self.file_format)))
            h_data_list.append([h_data, mc_cent])

        return h_data_list

    #-------------------------------------------------------------------------------------------
    # Plot ratio h1/h2
    def plot_ratio(self, h_pp, h_AA, outputFilename, xtitle, ytitle, eta_cut, cent, alpha_s, Q_switch,
                   label='Jet', R=None):

        # Create canvas
        cname = 'c_{}_{}_{}_{}_{}'.format(cent, alpha_s, Q_switch, eta_cut, R)
        c = ROOT.TCanvas(cname,cname,800,850)
        ROOT.SetOwnership(c, False) # For some reason this is necessary to avoid a segfault...some bug in ROOT or pyroot
                                    # Supposedly fixed in https://github.com/root-project/root/pull/3787
        c.cd()
        pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
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
        h_pp.SetMaximum(h_pp.GetMaximum()*10.)
        h_pp.SetMinimum(h_pp.GetMinimum()/10.)
            
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

        if 'Jet' in label:
            system3 = ROOT.TLatex(0.49 ,0.765, 'Anti-#it{k}_{T} #it{R} = 0.4 | #it{#eta}_{jet}| < ' + str(eta_cut))
        elif 'Hadron' in label:
            system3 = ROOT.TLatex(0.49 ,0.765, 'Charged particles | #it{#eta} | < ' + str(eta_cut))
            system3.SetNDC()
            system3.SetTextSize(0.044)
            system3.Draw()

        myLegend2pp = ROOT.TLegend(0.45,0.6,0.6,0.68)
        self.setupLegend(myLegend2pp,0.04)
        myLegend2pp.AddEntry(h_AA,'Pb-Pb {}%  (#alpha_{{s}}={}, Q_{{switch}}={})'.format(cent, alpha_s, Q_switch),'P')
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

        if self.debug_level > 0:
            c.SaveAs(outputFilename)
        
        return hRatio.Clone('{}_{}'.format(hRatio.GetName(), 'new'))
        
  
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
    print('Executing plot_results_PHYS_RAA.py...')
    print('')
    
    # Define arguments
    parser = argparse.ArgumentParser(description='Plot JETSCAPE events')
    parser.add_argument(
        '-c',
        '--configFile',
        action='store',
        type=str,
        metavar='configFile',
        default='config/PHYS_RAA.yaml',
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
