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
import pptx             # pip install python-pptx

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
        self.file_format = '.pdf'

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
        
        self.plot_hadron_correlation_observables(observable_type='hadron_correlations')
        
        self.plot_jet_observables(observable_type='inclusive_chjet')
        
        if 'inclusive_jet' in self.config:
            self.plot_jet_observables(observable_type='inclusive_jet')
            
        if 'semi_inclusive_chjet' in self.config:
            self.plot_semi_inclusive_chjet_observables(observable_type='semi_inclusive_chjet')
            
        if 'dijet' in self.config:
            self.plot_jet_observables(observable_type='dijet')
        
        # Generate pptx for convenience
        if self.file_format == '.png':
            self.generate_pptx()

    #-------------------------------------------------------------------------------------------
    # Plot hadron observables
    #-------------------------------------------------------------------------------------------
    def plot_hadron_observables(self, observable_type=''):
        print()
        print(f'Plot {observable_type} observables...')
                
        for observable, block in self.config[observable_type].items():
            for centrality_index,centrality in enumerate(block['centrality']):
        
                if 'hepdata' not in block:
                    continue
            
                # Initialize observable configuration
                self.suffix = ''
                self.init_observable(observable_type, observable, block, centrality)
                    
                # Plot observable
                self.plot_distribution_and_ratio(observable_type, observable, centrality)

    #-------------------------------------------------------------------------------------------
    # Plot hadron correlation observables
    #-------------------------------------------------------------------------------------------
    def plot_hadron_correlation_observables(self, observable_type=''):
        print()
        print(f'Plot {observable_type} observables...')

        for observable, block in self.config[observable_type].items():
            for centrality_index,centrality in enumerate(block['centrality']):

                if 'hepdata' not in block:
                    continue

                # Initialize observable configuration
                self.suffix = ''
                self.init_observable(observable_type, observable, block, centrality)
                  
                # Histogram observable
                self.plot_distribution_and_ratio(observable_type, observable, centrality)

    #-------------------------------------------------------------------------------------------
    # Histogram inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def plot_jet_observables(self, observable_type=''):
        print()
        print(f'Plot {observable_type} observables...')

        for observable, block in self.config[observable_type].items():
            for centrality_index,centrality in enumerate(block['centrality']):
                
                for self.jet_R in block['jet_R']:
                
                    # Optional: Loop through pt bins
                    for pt_bin in range(len(block['pt'])-1):

                        if len(block['pt']) > 2:
                            pt_suffix = f'_pt{pt_bin}'
                        else:
                            pt_suffix = ''
                            
                        # Optional: subobservable
                        subobservable_label_list = ['']
                        if 'kappa' in block:
                            subobservable_label_list = [f'_k{kappa}' for kappa in block['kappa']]
                        for subobservable_label in subobservable_label_list:
                        
                            # Set normalization
                            self_normalize = False
                            for x in ['mass', 'g', 'ptd', 'charge', 'mg', 'zg', 'tg']:
                                if x in observable:
                                    self_normalize = True
                        
                            if 'SoftDrop' in block:
                                for grooming_setting in block['SoftDrop']:
                                    if observable == 'tg_alice' and self.jet_R == 0.2 and grooming_setting['zcut'] == 0.4:
                                        continue
                                    else:
                                        print(f'      grooming_setting = {grooming_setting}')
                                        zcut = grooming_setting['zcut']
                                        beta = grooming_setting['beta']
                                        
                                        self.suffix = f'_R{self.jet_R}_zcut{zcut}_beta{beta}{subobservable_label}'
                                        if 'hepdata' not in block:
                                            continue
                            
                                        # Initialize observable configuration
                                        self.init_observable(observable_type, observable, block, centrality, pt_suffix=pt_suffix, self_normalize=self_normalize)
                                  
                                        # Plot observable
                                        self.plot_distribution_and_ratio(observable_type, observable, centrality, pt_suffix)
                                
                            else:

                                self.suffix = f'_R{self.jet_R}{subobservable_label}'
                                if 'hepdata' not in block:
                                    continue

                                # Initialize observable configuration
                                self.init_observable(observable_type, observable, block, centrality, pt_suffix=pt_suffix, self_normalize=self_normalize)
                          
                                # Plot observable
                                self.plot_distribution_and_ratio(observable_type, observable, centrality, pt_suffix)

    #-------------------------------------------------------------------------------------------
    # Histogram semi-inclusive jet observables
    #-------------------------------------------------------------------------------------------
    def plot_semi_inclusive_chjet_observables(self, observable_type=''):
        print()
        print(f'Plot {observable_type} observables...')
        
        for observable, block in self.config[observable_type].items():
            for centrality_index,centrality in enumerate(block['centrality']):
                    
                for self.jet_R in block['jet_R']:
                    
                    self.suffix = f'_R{self.jet_R}'
                    if 'hepdata' not in block:
                        continue
                        
                    # Set normalization
                    self_normalize = False
                    for x in ['nsubjettiness']:
                        if x in observable:
                            self_normalize = True

                    # Initialize observable configuration
                    self.init_observable(observable_type, observable, block, centrality, self_normalize=self_normalize)
              
                    # Plot observable
                    self.plot_distribution_and_ratio(observable_type, observable, centrality)

    #-------------------------------------------------------------------------------------------
    # Initialize a single observable's config
    #-------------------------------------------------------------------------------------------
    def init_observable(self, observable_type, observable, block, centrality, pt_suffix='', self_normalize=False):
    
        # Initialize an empty dict containing relevant info
        self.observable_settings = {}
                            
        # Common settings
        self.xtitle = block['xtitle']
        if 'eta_cut' in block:
            self.eta_cut = block['eta_cut']
        if 'pt' in block:
            self.pt = block['pt']
        if 'eta_R' in block:
            self.eta_R = block['eta_R']
            self.eta_cut = np.round(self.eta_R - self.jetR, decimals=1)
        if 'c_ref' in block:
            index = block['jet_R'].index(self.jet_R)
            self.c_ref = block['c_ref'][index]
        if 'low_trigger_range' in block:
            low_trigger_range = block['low_trigger_range']
        if 'high_trigger_range' in block:
            high_trigger_range = block['high_trigger_range']
        if 'trigger_range' in block:
            trigger_range = block['trigger_range']
        if 'logy' in block:
            self.logy = block['logy']
        else:
            self.logy = False
        
        if 'ytitle_pp' in block:
            self.ytitle = block['ytitle_pp']
        else:
            self.ytitle = ''
        if 'y_min_pp' in block:
            self.y_min = float(block['y_min_pp'])
            self.y_max = float(block['y_max_pp'])
        else:
            self.y_min = 0.
            self.y_max = 1.
        if 'y_ratio_min' in block:
            self.y_ratio_min = block['y_ratio_min']
            self.y_ratio_max = block['y_ratio_max']
        else:
            self.y_ratio_min = 0.
            self.y_ratio_max = 1.99
        if 'skip_pp_ratio' in block:
            self.skip_pp_ratio = block['skip_pp_ratio']
        else:
            self.skip_pp_ratio = False
        if 'scale_by' in block:
            self.scale_by = block['scale_by']
        else:
            self.scale_by = None
        
        # Initialize data
        if f'hepdata' in block:
            self.observable_settings['data_distribution'] = self.plot_utils.tgraph_from_hepdata(block, self.sqrts, observable_type, observable, suffix=self.suffix, pt_suffix=pt_suffix)
        else:
            self.observable_settings['data_distribution'] = None

        # Initialize JETSCAPE
        keys = [key.ReadObj().GetTitle() for key in self.input_file.GetListOfKeys()]
        if 'semi_inclusive' not in observable_type:
            self.hname = f'h_{observable_type}_{observable}{self.suffix}_{centrality}{pt_suffix}'
            if self.hname in keys:
                h_jetscape = self.input_file.Get(self.hname)
                h_jetscape.SetDirectory(0)
            else:
                h_jetscape = None
            self.observable_settings['jetscape_distribution'] = h_jetscape
        else:
            if self.sqrts == 2760: # Delta recoil
                hname_low_trigger = f'h_{observable_type}_{observable}_R{self.jet_R}_lowTrigger_{centrality}'
                hname_high_trigger = f'h_{observable_type}_{observable}_R{self.jet_R}_highTrigger_{centrality}'
                hname_ntrigger = f'h_{observable_type}_alice_trigger_pt{observable}_{centrality}'
                self.hname = f'h_{observable_type}_{observable}_R{self.jet_R}_{centrality}'
                if hname_low_trigger in keys and hname_high_trigger in keys and hname_ntrigger in keys:
                    h_jetscape_low = self.input_file.Get(hname_low_trigger)
                    h_jetscape_low.SetDirectory(0)
                    h_jetscape_high = self.input_file.Get(hname_high_trigger)
                    h_jetscape_high.SetDirectory(0)
                    h_jetscape_ntrigger = self.input_file.Get(hname_ntrigger)
                    h_jetscape_ntrigger.SetDirectory(0)
                    
                    low_trigger = (low_trigger_range[0]+low_trigger_range[1])/2
                    high_trigger = (high_trigger_range[0]+high_trigger_range[1])/2
                    n_trig_high = h_jetscape_ntrigger.GetBinContent(h_jetscape_ntrigger.FindBin(high_trigger))
                    n_trig_low = h_jetscape_ntrigger.GetBinContent(h_jetscape_ntrigger.FindBin(low_trigger))
                    h_jetscape_high.Scale(1./n_trig_high)
                    h_jetscape_low.Scale(1./n_trig_low)
                    self.observable_settings['jetscape_distribution'] = h_jetscape_high.Clone('h_delta_recoil')
                    self.observable_settings['jetscape_distribution'].Add(h_jetscape_low, -1)
                else:
                    self.observable_settings['jetscape_distribution'] = None
            elif self.sqrts == 200:
                hname = f'{observable_type}_{observable}_R{jet_R}'
                hname_ntrigger = f'{observable_type}_star_trigger_pt'
                if hname in keys and hname_trigger in keys:
                    self.observable_settings['jetscape_distribution'] = self.input_file.Get(hname)
                    self.observable_settings['jetscape_distribution'].SetDirectory(0)
                    h_jetscape_ntrigger = self.input_file.Get(hname_n_trigger)
                    h_jetscape_ntrigger.SetDirectory(0)
                    
                    trigger = (trigger_range[0]+trigger_range[1])/2
                    n_trig = h_jetscape_ntrigger.GetBinContent(h_jetscape_ntrigger.FindBin(trigger))
                    self.observable_settings['jetscape_distribution'].Scale(1./n_trig)
                else:
                    self.observable_settings['jetscape_distribution'] = None
        
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
        # Note: If we divide by n_events, then JETSCAPE distribution gives cross-section (in b)
        if self.observable_settings['jetscape_distribution']:
            n_events = self.input_file.Get('h_n_events').GetBinContent(1)
            self.observable_settings['jetscape_distribution'].Scale(1./n_events)
            if self.sqrts == 2760:
                if observable_type == 'hadron':
                    if observable == 'pt_ch_alice':
                        sigma_inel = 0.0618
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1./sigma_inel)
                    elif observable == 'pt_pi_alice':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                    elif observable == 'pt_pi0_alice':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                    elif observable == 'pt_ch_atlas':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*np.pi))
                        self.observable_settings['jetscape_distribution'].Scale(1.e3) # convert to mb
                    elif observable == 'pt_ch_cms':
                        L_int = 230. # (in nb)
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*np.pi))
                        self.observable_settings['jetscape_distribution'].Scale(L_int)
                elif observable_type == 'inclusive_jet':
                    if observable == 'pt_alice':
                        sigma_inel = 0.0621
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1./sigma_inel)
                    if observable == 'pt_atlas':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1.e9) # convert to nb
                    if observable == 'pt_cms':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                        self.observable_settings['jetscape_distribution'].Scale(1.e9) # convert to nb
                    if observable in ['Dz_atlas', 'Dpt_atlas']:
                        hname = f'h_{observable_type}_{observable}{self.suffix}_Njets_{centrality}{pt_suffix}'
                        h_njets = self.input_file.Get(self.hname)
                        h_njets.SetDirectory(0)
                        n_jets = h_njets.GetBinContent(1)
                        if n_jets > 0.:
                            self.observable_settings['jetscape_distribution'].Scale(1.*n_events/n_jets)
                        else:
                            print('WARNING: N_jets = 0')
                elif observable_type == 'inclusive_chjet':
                    if observable == 'pt_alice':
                        self.observable_settings['jetscape_distribution'].Scale(1./(2*self.eta_cut))
                elif observable_type == 'semi_inclusive_chjet':
                    if observable in ['IAA_alice', 'dphi_alice']:
                        self.observable_settings['jetscape_distribution'].Scale(n_events)

            # Scale by bin-dependent factor (approximation for pp QA)
            if self.scale_by:
                nBins = self.observable_settings['jetscape_distribution'].GetNbinsX()
                for bin in range(1, nBins+1):
                    h_x = self.observable_settings['jetscape_distribution'].GetBinCenter(bin)
                    h_y = self.observable_settings['jetscape_distribution'].GetBinContent(bin)
                    h_error = self.observable_settings['jetscape_distribution'].GetBinError(bin)
                    if self.scale_by == '1/pt':
                        scale_factor = 1./h_x
                    self.observable_settings['jetscape_distribution'].SetBinContent(bin, h_y*scale_factor)
                    self.observable_settings['jetscape_distribution'].SetBinError(bin, h_error*scale_factor)

            # Bin width
            self.observable_settings['jetscape_distribution'].Scale(1., 'width')
            
            # Self-normalization
            if self_normalize:
                if observable in ['zg', 'theta_g']:
                    min_bin = 0
                else:
                    min_bin = 1
                nbins = self.observable_settings['jetscape_distribution'].GetNbinsX()
                integral = self.observable_settings['jetscape_distribution'].Integral(min_bin, nbins, 'width')
                if integral > 0.:
                    self.observable_settings['jetscape_distribution'].Scale(1./integral)
                else:
                    print('WARNING: integral for self-normalization is 0')
                    
            if observable == 'nsubjettiness_alice':
                print(self.observable_settings['jetscape_distribution'].Integral())
                
        # Form ratio of JETSCAPE to data
        if self.observable_settings['data_distribution'] and self.observable_settings['jetscape_distribution'] and not self.observable_settings['jetscape_distribution'].InheritsFrom(ROOT.TH2.Class()) and not self.skip_pp_ratio:
            self.observable_settings['ratio'] = self.plot_utils.divide_histogram_by_tgraph(self.observable_settings['jetscape_distribution'],
                                                                              self.observable_settings['data_distribution'])
                                                                              
        else:
            self.observable_settings['ratio'] = None

    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, observable_type, observable, centrality, pt_suffix='', logy = False):
    
        if not self.observable_settings['jetscape_distribution']:
            return
            
        self.observable_settings['jetscape_distribution']
    
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
        if self.logy:
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
        if self.observable_settings['data_distribution']:
            self.observable_settings['data_distribution'].SetMarkerSize(self.marker_size)
            self.observable_settings['data_distribution'].SetMarkerStyle(self.data_marker)
            self.observable_settings['data_distribution'].SetMarkerColor(self.data_color)
            self.observable_settings['data_distribution'].SetLineStyle(self.line_style)
            self.observable_settings['data_distribution'].SetLineWidth(self.line_width)
            self.observable_settings['data_distribution'].SetLineColor(self.data_color)
            self.observable_settings['data_distribution'].Draw('PE Z same')
            legend.AddEntry(self.observable_settings['data_distribution'], 'Data', 'PE')
        
        legend.Draw()
        
        # Draw ratio
        pad2.cd()
        if self.observable_settings['ratio']:
            self.observable_settings['ratio'].SetFillColor(self.jetscape_color)
            self.observable_settings['ratio'].SetFillColorAlpha(self.jetscape_color, self.alpha)
            self.observable_settings['ratio'].SetFillStyle(1001)
            self.observable_settings['ratio'].SetMarkerSize(0.)
            self.observable_settings['ratio'].SetMarkerStyle(0)
            self.observable_settings['ratio'].SetLineWidth(0)
            self.observable_settings['ratio'].Draw('E3 same')
        
        # Draw data uncertainties
        if self.observable_settings['data_distribution']:
            data_ratio = self.plot_utils.divide_tgraph_by_tgraph(self.observable_settings['data_distribution'],
                                                                 self.observable_settings['data_distribution'])
            data_ratio.Draw('PE Z same')

        line = ROOT.TLine(self.bins[0], 1, self.bins[-1], 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        
        pad1.cd()
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        x = 0.25
        text_latex.SetTextSize(0.065)
        text = f'#bf{{{observable_type}_{observable}}} #sqrt{{#it{{s}}}} = {self.sqrts/1000.} TeV'
        text_latex.DrawLatex(x, 0.83, text)
        text = f'{centrality} {self.suffix} {pt_suffix}'
        text_latex.DrawLatex(x, 0.73, text)

        c.SaveAs(os.path.join(self.output_dir, f'{self.hname}{self.file_format}'))
        c.Close()

    #-------------------------------------------------------------------------------------------
    # Generate pptx of one plot per slide, for convenience
    #-------------------------------------------------------------------------------------------
    def generate_pptx(self):
    
        # Create a blank presentation
        p = pptx.Presentation()
        
        # Set slide layouts
        title_slide_layout = p.slide_layouts[0]
        blank_slide_layout = p.slide_layouts[6]
        
        # Make a title slide
        slide = p.slides.add_slide(title_slide_layout)
        title = slide.shapes.title
        title.text = f'QA for {self.sqrts/1000.} TeV'
        author = slide.placeholders[1]
        author.text = 'STAT WG'
        
        # Loop through all output files and plot
        files = [f for f in os.listdir(self.output_dir) if f.endswith(self.file_format)]
        for file in sorted(files):
            img = os.path.join(self.output_dir, file)
            slide = p.slides.add_slide(blank_slide_layout)
            slide.shapes.add_picture(img, left=pptx.util.Inches(2.),
                                          top=pptx.util.Inches(1.),
                                          width=pptx.util.Inches(5.))
            
        p.save(os.path.join(self.output_dir, 'results.pptx'))

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
