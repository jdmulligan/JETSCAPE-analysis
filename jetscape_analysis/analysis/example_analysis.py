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

from jetscape_analysis.analysis.event import event_hepmc
from jetscape_analysis.base import common_base

################################################################
class example_analysis(common_base.common_base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file="", input_file="", output_dir="", bin="", **kwargs):
        super(example_analysis, self).__init__(**kwargs)
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
        with open(self.config_file, "r") as stream:
            config = yaml.safe_load(stream)

        self.debug_level = config["debug_level"]
        self.jetR_list = config["jetR"]
        self.min_jet_pt = config["min_jet_pt"]
        self.abs_jet_eta_max = config["abs_jet_eta_max"]
        self.event_id = 0

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_output_objects(self):

        # Event histograms
        self.hNevents = ROOT.TH1F("hNevents", "hNevents", 20, 0, 20)
        self.hCrossSection = ROOT.TH1F("hCrossSection", "hCrossSection", 20, 0, 20)

        # Hadron histograms
        self.hHadronN = ROOT.TH1F("hHadronN", "hHadronN", 1000, 0, 1000)
        self.hHadronPt = ROOT.TH1F("hHadronPt", "hHadronPt", 300, 0.0, 300.0)
        self.hHadronPID = ROOT.TH1F("hHadronPID", "hHadronPID", 10000, -5000, 5000)
        self.hHadronEtaPhi = ROOT.TH2F("hHadronEtaPhi", "hHadronEtaPhi", 100, -5, 5, 100, -3.2, 3.2)

        # Final-state parton histograms
        self.hPartonN = ROOT.TH1F("hPartonN", "hPartonN", 1000, 0, 1000)
        self.hPartonPt = ROOT.TH1F("hPartonPt", "hPartonPt", 300, 0.0, 300.0)
        self.hPartonPID = ROOT.TH1F("hPartonPID", "hPartonPID", 10000, -5000, 5000)
        self.hPartonEtaPhi = ROOT.TH2F("hPartonEtaPhi", "hPartonEtaPhi", 100, -5, 5, 100, -3.2, 3.2)

        # Jet histograms
        for jetR in self.jetR_list:

            name = "hJetPt_R{}".format(jetR)
            h = ROOT.TH1F(name, name, 300, 0, 300)
            setattr(self, name, h)

            name = "hJetEtaPhi_R{}".format(jetR)
            h = ROOT.TH2F(name, name, 100, -5, 5, 100, -6.28, 6.28)
            setattr(self, name, h)

    # ---------------------------------------------------------------
    # Analyze a single event
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Print and store basic event info
        self.get_event_info(event)

        # Get list of hadrons from the event, and fill some histograms
        hadrons = event.hadrons()
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
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(self.abs_jet_eta_max - jetR)
            if self.debug_level > 0:
                print("jet definition is:", jet_def)
                print("jet selector is:", jet_selector, "\n")

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

        # Get cross-section
        # cross_section = event.GenCrossSection
        # fCrossSection = cross_section.xsec(0)
        # xsec_error = cross_section.xsec_err(0)

        # Print some basic info for first event
        # if self.event_id == 0:

        # print('xsec: {} +/- {} pb'.format(fCrossSection, xsec_error))

        # Get heavy ion attributes
        # heavy_ion = event.GenHeavyIon
        # nColl = heavy_ion.Ncoll
        # nPart = heavy_ion.Npart_proj
        # eventPlaneAngle = heavy_ion.event_plane_angle
        # print('NColl = {}, NPart = {}, EP-angle = {}'.format(nColl, nPart, eventPlaneAngle))

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

            if abs(eta) < self.abs_jet_eta_max:
                self.hHadronPt.Fill(pt)
                self.hHadronPID.Fill(pid)
                self.hHadronEtaPhi.Fill(eta, phi)

        self.hHadronN.Fill(len(hadrons))

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

            if abs(eta) < self.abs_jet_eta_max:
                self.hPartonPt.Fill(pt)
                self.hPartonPID.Fill(pid)
                self.hPartonEtaPhi.Fill(eta, phi)

        self.hPartonN.Fill(len(partons))

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

            h = getattr(self, "hJetPt_R{}".format(jetR))
            h.Fill(jet_pt)

            h = getattr(self, "hJetEtaPhi_R{}".format(jetR))
            h.Fill(jet_eta, jet_phi)

    # ---------------------------------------------------------------
    # Save all ROOT histograms and trees to file
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Fill cross-section with last event's value, which is most accurate
        # self.hCrossSection(fCrossSection/(1e9))
        self.hNevents.SetBinContent(self.bin + 1, self.event_id)

        # Save output objects
        outputfilename = os.path.join(self.output_dir, "AnalysisResults.root")
        fout = ROOT.TFile(outputfilename, "recreate")
        fout.cd()
        for attr in dir(self):

            obj = getattr(self, attr)

            # Write all ROOT histograms and trees to file
            types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
            if isinstance(obj, types):
                obj.Write()

        fout.Close()

    # ---------------------------------------------------------------
    # Remove periods from a label
    # ---------------------------------------------------------------
    def remove_periods(self, text):

        string = str(text)
        return string.replace(".", "")
