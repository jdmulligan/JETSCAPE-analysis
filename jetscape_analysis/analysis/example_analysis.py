#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE output file

  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# General
import os

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjext
import ROOT
import yaml
from array import *

from jetscape_analysis.analysis.event import event_hepmc
from jetscape_analysis.base import common_base

################################################################
class ExampleAnalysis(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', bin='', **kwargs):
        super(ExampleAnalysis, self).__init__(**kwargs)
        self.config_file = config_file
        self.input_file = input_file
        self.output_dir = output_dir
        self.bin = bin

        self.initialize_config()
        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.parameter_scan_dict = config['parameter_scan']
        self.n_pt_hat_bins = len(self.parameter_scan_dict['pt_hat_bins']['values']) - 1

        self.debug_level = config['debug_level']
        self.reader = config['reader']
        self.event_id = 0
        
        self.min_track_pt = config['min_track_pt']
        self.abs_track_eta_max = config['abs_track_eta_max']
        
        self.jetR_list = config['jetR']
        self.min_jet_pt = config['min_jet_pt']

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_output_objects(self):

        # Event histograms
        self.hNevents = ROOT.TH1F('hNevents', 'hNevents', self.n_pt_hat_bins, 0, self.n_pt_hat_bins)
        self.hCrossSection = ROOT.TH1F('hCrossSection', 'hCrossSection', self.n_pt_hat_bins, 0, self.n_pt_hat_bins)

        # Hadron histograms
        hname = 'hChHadronPt_eta'
        pt_bins = [9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 140.0, 165.0, 250.0, 400.0]
        n_pt_bins = len(pt_bins) - 1
        pt_bin_array = array('d', pt_bins)
        eta_bins = [-5., -3., -1., 1., 3., 5.]
        n_eta_bins = len(eta_bins) - 1
        eta_bin_array = array('d', eta_bins)
        setattr(self, hname, ROOT.TH2F(hname, hname, n_pt_bins, pt_bin_array, n_eta_bins, eta_bin_array))
        
        hname = 'hD0Pt_eta'
        pt_bins = [2., 3., 4., 5., 6., 8., 10., 12.5, 15., 20., 25., 30., 40., 60., 100.]
        n_pt_bins = len(pt_bins) - 1
        pt_bin_array = array('d', pt_bins)
        eta_bins = [-5., -3., -1., 1., 3., 5.]
        n_eta_bins = len(eta_bins) - 1
        eta_bin_array = array('d', eta_bins)
        setattr(self, hname, ROOT.TH2F(hname, hname, n_pt_bins, pt_bin_array, n_eta_bins, eta_bin_array))
        
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
    # Analyze a single event
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Print and store basic event info
        self.get_event_info(event)

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
    # Get event info
    # ---------------------------------------------------------------
    def get_event_info(self, event):

        # Print some basic info for first event
        #if self.event_id == 0:

            # Get heavy ion attributes
            #heavy_ion = event.heavy_ion()
            # However it seems that pyhepmc_ng doesn't implement any of these...
            #print(dir(heavy_ion))
            #nColl = heavy_ion.Ncoll
            #nPart = heavy_ion.Npart_proj
            #eventPlaneAngle = heavy_ion.event_plane_angle
            #print('NColl = {}, NPart = {}, EP-angle = {}'.format(nColl, nPart, eventPlaneAngle))

        self.event_id += 1

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

            #  Fill charged particle histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-, e+, mu+)
            #  (assuming weak strange decays are off, but charm decays are on)
            if abs(pid)==211 or abs(pid)==321 or abs(pid)==2212 or abs(pid)==3222 or abs(pid)==3112 or abs(pid)==3312 or abs(pid)==3334 or abs(pid)==11 or abs(pid)==13:
                getattr(self, 'hChHadronPt_eta').Fill(pt, eta)
                getattr(self, 'hChHadronEtaPhi').Fill(eta, phi)

            if abs(eta) < self.abs_track_eta_max:
            
                hname = 'hHadronPID_eta{}'.format(self.abs_track_eta_max)
                getattr(self, hname).Fill(pid)
                
            if abs(pid) == 421:
            
                getattr(self, 'hD0Pt_eta').Fill(pt, eta)

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
    # Fill hadrons into vector of fastjet pseudojets
    # ---------------------------------------------------------------
    def fill_fastjet_constituents(self, hadrons):

        px = [hadron.momentum.px for hadron in hadrons]
        py = [hadron.momentum.py for hadron in hadrons]
        pz = [hadron.momentum.pz for hadron in hadrons]
        e = [hadron.momentum.e for hadron in hadrons]
        
        # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of px,py,pz,e
        fj_particles = fjext.vectorize_px_py_pz_e(px, py, pz, e)
        
        return fj_particles
    
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

    # ---------------------------------------------------------------
    # Save all ROOT histograms and trees to file
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Fill cross-section with last event's value, which is most accurate
        xsec = self.cross_section()
        self.hCrossSection.SetBinContent(self.bin+1, xsec)
        
        # Set N events
        self.hNevents.SetBinContent(self.bin + 1, self.event_id)

        # Save output objects
        outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
        fout = ROOT.TFile(outputfilename, 'recreate')
        fout.cd()
        for attr in dir(self):

            obj = getattr(self, attr)

            # Write all ROOT histograms and trees to file
            types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
            if isinstance(obj, types):
                obj.Write()

        fout.Close()

    # ---------------------------------------------------------------
    # Get cross-section from last event JETSCAPE output file
    #
    # It seems that pyhepmc_ng doesn't contain GenCrossSection, so we need to find it manually
    # Similarly, we find the cross-section manually for ascii format.
    # ---------------------------------------------------------------
    def cross_section(self):
        
        # Fill array of cross-sections
        cross_sections = []
        with open(self.input_file, 'r') as infile:
            for line in infile:
            
                if self.reader == 'hepmc':
                    if 'GenCrossSection' in line:
                        split = line.split()
                        xsec = float(split[3]) / 1e9
                        cross_sections.append(xsec)
            
                elif self.reader == 'ascii':
                    if 'sigmaGen' in line:
                        split = line.split()
                        xsec = float(split[2])
                        cross_sections.append(line)
                        
        # Return cross-section with last event's value, which is most accurate
        return cross_sections[-1]

    # ---------------------------------------------------------------
    # Remove periods from a label
    # ---------------------------------------------------------------
    def remove_periods(self, text):

        string = str(text)
        return string.replace('.', '')
