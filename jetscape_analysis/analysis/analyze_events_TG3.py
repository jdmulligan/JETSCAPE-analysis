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

        #------------------------------------------------------
        # Charged particle parameters
        self.charged_particle_eta_cut = config['charged_particle_eta_cut']

        self.file_CMS_hadron_0_5 = config['CMS_hadron_0_5']
        self.file_CMS_hadron_5_10 = config['CMS_hadron_5_10']
        self.file_CMS_hadron_30_50 = config['CMS_hadron_30_50']
        self.file_ATLAS_hadron = config['ATLAS_hadron']
        self.file_ALICE_hadron = config['ALICE_hadron']

        # Get binnings from data
        f = ROOT.TFile(self.file_CMS_hadron_0_5, 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.bins_CMS_hadron = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.file_ALICE_hadron, 'READ')
        dir = f.Get('Table 8')
        h = dir.Get('Hist1D_y1')
        self.bins_ALICE_hadron = np.array(h.GetXaxis().GetXbins())
        f.Close()

        self.bins_ATLAS_hadron = np.array([5.5, 7, 8, 9, 10, 11, 13, 15, 17, 20, 22, 25, 30, 38, 48, 59, 74, 90, 120, 150, 200, 250, 300])

        #------------------------------------------------------
        # Jet parameters
        self.use_charged_jets = config["use_charged_jets"]
        self.min_track_pt = config['min_track_pt']
        self.jetR_list = config['jetR']
        self.min_jet_pt = config['min_jet_pt']
        self.jet_eta_cut_04 = config['jet_eta_cut_04']
        self.jet_eta_cut_02 = config['jet_eta_cut_02']

        self.file_CMS_jet = config['CMS_jet']
        self.file_ATLAS_jet_0_10 = config['ATLAS_jet_0_10']
        self.file_ATLAS_jet_30_40 = config['ATLAS_jet_30_40']
        self.file_ATLAS_jet_40_50 = config['ATLAS_jet_40_50']
        self.file_ALICE_jet_0_10_R02 = config['ALICE_jet_0_10_R02']
        self.file_ALICE_jet_0_10_R04 = config['ALICE_jet_0_10_R04']

        # Get binnings from data
        self.bins_CMS_jet = np.array([250., 300., 400., 500., 1000.])

        f = ROOT.TFile(self.file_ATLAS_jet_0_10, 'READ')
        dir = f.Get('Table 19')
        h = dir.Get('Hist1D_y1')
        self.bins_ATLAS_jet_0_10 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.file_ATLAS_jet_30_40, 'READ')
        dir = f.Get('Table 22')
        h = dir.Get('Hist1D_y1')
        self.bins_ATLAS_jet_30_40 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.file_ATLAS_jet_40_50, 'READ')
        dir = f.Get('Table 23')
        h = dir.Get('Hist1D_y1')
        self.bins_ATLAS_jet_40_50 = np.array(h.GetXaxis().GetXbins())
        f.Close()

        f = ROOT.TFile(self.file_ALICE_jet_0_10_R02, 'READ')
        dir = f.Get('Table 30')
        h = dir.Get('Hist1D_y1')
        self.bins_ALICE_jet = np.array(h.GetXaxis().GetXbins())
        f.Close()

        # Angularity
        # TODO: Need the HEPdata.
        self.bins_ALICE_angularity = ...

        # Jet mass
        # TODO: Need the HEPdata.
        self.bins_ALICE_jet_mass = ...

        # Substructure
        self.soft_drop_zcut = config["soft_drop_zcut"]
        self.soft_drop_beta = config["soft_drop_beta"]
        # TODO: Need binning
        self.bins_ALICE_soft_drop_zg = ...
        self.bins_ALICE_soft_drop_theta_g = ...

        # Hadron-jet
        self.hjet_low_trigger_range = config["hjet_low_trigger_range"]
        self.hjet_high_trigger_range = config["hjet_high_trigger_range"]
        # TODO: Need the HEPdata.
        self.bins_ALICE_hjet_yield = ...
        self.bins_ALICE_hjet_delta_phi = ...


    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):

        # Hadron histograms
        hname = 'hChargedPt_CMS'
        h = ROOT.TH1F(hname, hname, len(self.bins_CMS_hadron)-1, self.bins_CMS_hadron)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hChargedPt_ATLAS'
        h = ROOT.TH1F(hname, hname, len(self.bins_ATLAS_hadron)-1, self.bins_ATLAS_hadron)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hChargedPt_ALICE'
        h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_hadron)-1, self.bins_ALICE_hadron)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hChargedPt_Recoils'
        h = ROOT.TH1F(hname, hname, 3000, 0, 300)
        h.Sumw2()
        setattr(self, hname, h)

        # Jet histograms
        for jetR in self.jetR_list:

            hname = 'hJetPt_CMS_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, len(self.bins_CMS_jet)-1, self.bins_CMS_jet)
            h.Sumw2()
            setattr(self, hname, h)

            hname = 'hJetPt_ALICE_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_jet)-1, self.bins_ALICE_jet)
            h.Sumw2()
            setattr(self, hname, h)

            hname = 'hJetPt_ALICE_no_ptlead_cut_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_jet)-1, self.bins_ALICE_jet)
            h.Sumw2()
            setattr(self, hname, h)

            # Angularity
            hname = f"hAngularity_R{jetR}"
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_angularity)-1, self.bins_ALICE_angularity)
            h.Sumw2()
            setattr(self, hname, h)

            # Jet mass
            hname = f"hJetMass_R{jetR}"
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_jet_mass)-1, self.bins_ALICE_jet_mass)
            h.Sumw2()
            setattr(self, hname, h)

            # Substructure
            # zg
            hname = f"hSoftDrop_zg_R{jetR}"
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_soft_drop_zg)-1, self.bins_ALICE_soft_drop_zg)
            h.Sumw2()
            setattr(self, hname, h)
            # theta_g
            hname = f"hSoftDrop_theta_g_R{jetR}"
            h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_soft_drop_theta_g)-1, self.bins_ALICE_soft_drop_theta_g)
            h.Sumw2()
            setattr(self, hname, h)

            # h-jet
            for hist_label in ["low", "high"]:
                # Yield
                hname = f"hRecoilJetYield_{hist_label}Trigger"
                h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_hjet_yield)-1, self.bins_ALICE_hjet_yield)
                h.Sumw2()
                setattr(self, hname, h)
                # Delta Phi
                hname = f"hRecoilJetDeltaPhi_{hist_label}Trigger"
                h = ROOT.TH1F(hname, hname, len(self.bins_ALICE_hjet_delta_phi)-1, self.bins_ALICE_hjet_delta_phi)
                h.Sumw2()
                setattr(self, hname, h)

            hname = 'hJetPt_recoils_R{}'.format(jetR)
            h = ROOT.TH2F(hname, hname, 100, 0, 1000, 300, 0, 300)
            h.GetXaxis().SetTitle('jet pt')
            h.GetYaxis().SetTitle('recoil pt')
            h.Sumw2()
            setattr(self, hname, h)

        hname = 'hJetPt_ATLAS_binning0_R{}'.format(0.4)
        h = ROOT.TH1F(hname, hname, len(self.bins_ATLAS_jet_0_10)-1, self.bins_ATLAS_jet_0_10)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hJetPt_ATLAS_binning1_R{}'.format(0.4)
        h = ROOT.TH1F(hname, hname, len(self.bins_ATLAS_jet_30_40)-1, self.bins_ATLAS_jet_30_40)
        h.Sumw2()
        setattr(self, hname, h)

        hname = 'hJetPt_ATLAS_binning2_R{}'.format(0.4)
        h = ROOT.TH1F(hname, hname, len(self.bins_ATLAS_jet_40_50)-1, self.bins_ATLAS_jet_40_50)
        h.Sumw2()
        setattr(self, hname, h)

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):
        # Only used charged jets.
        if self.use_charged_jets:
            # NOTE: This is super inefficient - we can seriously accelerate this with numba.
            #       But this apparently works for now, so we leave it as is for now.
            # Create an all false mask. We'll fill it in with the charged constituents
            mask = np.ones(len(event)) < 0
            for i, particle in enumerate(event):
                # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if abs(particle.particle_ID) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    mask[i] = True
            # Can't apply directly to the event because not everything has the same depth
            # (eg. event plane angle is attached). So we extract out the columns where we
            # can do so, and pass that to fastjet. In that case, we can apply the mask to
            # the entire array.
            fj_event = event[["px", "py", "pz", "E", "particle_ID", "status"]][mask]
        else:
            fj_event = event

        # Create list of fastjet::PseudoJets (separately for jet shower particles and holes)
        fj_hadrons_positive = self.fill_fastjet_constituents(fj_event, select_status='+')
        fj_hadrons_negative = self.fill_fastjet_constituents(fj_event, select_status='-')

        # Fill hadron histograms for jet shower particles
        self.fill_hadron_histograms(fj_hadrons_positive, status='+')
        self.fill_hadron_histograms(fj_hadrons_negative, status='-')

        # Loop through specified jet R
        for jetR in self.jetR_list:

            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(5.)
            if self.debug_level > 0:
                print('jet definition is:', jet_def)
                print('jet selector is:', jet_selector, '\n')

            # Do jet finding
            cs = fj.ClusterSequence(fj_hadrons_positive, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # Fill some jet histograms
            self.fill_jet_histograms(jets_selected, fj_hadrons_positive, fj_hadrons_negative, jetR)

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
                getattr(self, 'hChargedPt_Recoils').Fill(pt)
                continue

            # CMS
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.charged_particle_eta_cut[0]:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_CMS').Fill(pt)

            # ATLAS
            # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            # Exclude e+, mu+ (11, 13)
            if abs(eta) < self.charged_particle_eta_cut[1]:
                if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_ATLAS').Fill(pt)

            # ALICE
            # Fill charged hadron histograms (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(eta) < self.charged_particle_eta_cut[2]:
                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    getattr(self, 'hChargedPt_ALICE').Fill(pt)

    # ---------------------------------------------------------------
    # Fill jet histograms
    # ---------------------------------------------------------------
    def fill_jet_histograms(self, jets, fj_hadrons_positive, fj_hadrons_negative, jetR):

        for jet in jets:

            # Get jet pt (from shower hadrons only)
            jet_pt_uncorrected = jet.pt()

            # Sum the negative recoils within R
            holes_in_jet = []
            negative_pt = 0.
            for hadron in fj_hadrons_negative:
                if jet.delta_R(hadron) < jetR:
                    negative_pt += hadron.pt()
                    holes_in_jet.append(hadron)

            # Compute corrected jet pt, and fill histograms
            jet_pt = jet_pt_uncorrected - negative_pt

            # Select eta cut
            if jetR == 0.2:
                jet_eta_cut = self.jet_eta_cut_02
            elif jetR == 0.4:
                jet_eta_cut = self.jet_eta_cut_04

            # CMS
            if abs(jet.eta()) < jet_eta_cut[0]:
                getattr(self, 'hJetPt_CMS_R{}'.format(jetR)).Fill(jet_pt)

            # ATLAS
            if jetR == 0.4:
                if abs(jet.rap()) < jet_eta_cut[1]:
                    getattr(self, 'hJetPt_ATLAS_binning0_R{}'.format(jetR)).Fill(jet_pt)
                    getattr(self, 'hJetPt_ATLAS_binning1_R{}'.format(jetR)).Fill(jet_pt)
                    getattr(self, 'hJetPt_ATLAS_binning2_R{}'.format(jetR)).Fill(jet_pt)

            # ALICE
            if abs(jet.eta()) < jet_eta_cut[2]:

                # Check leading track requirement
                if jetR == 0.2:
                    min_leading_track_pt = 5.
                elif jetR == 0.4:
                    min_leading_track_pt = 7.

                accept_jet = False
                for constituent in jet.constituents():
                    if constituent.pt() > min_leading_track_pt:
                        # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                        if abs(constituent.user_index()) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                            accept_jet = True

                if accept_jet:
                    getattr(self, 'hJetPt_ALICE_R{}'.format(jetR)).Fill(jet_pt)
                getattr(self, 'hJetPt_ALICE_no_ptlead_cut_R{}'.format(jetR)).Fill(jet_pt)

                # Angularity
                # No LTB, so take all jets within acceptance.
                # NOTE: Charged jet pt
                if 40 < jet_pt < 60:
                    value = 0
                    # Intentionally use uncorrected jet pt here because our approach is to
                    # calculate the unmeasured angularity, and then separately calculate a
                    # value from the holes to be subtracted.
                    for constituent in jet.constituents():
                        value += constituent.pt() / jet.pt() * constituent.delta_R(jet)
                    # Find holes within the jet cone, and subtract their contribution
                    for hadron in holes_in_jet:
                        value -= hadron.pt() / jet.pt() * hadron.delta_R(jet)
                    getattr(self, f"hAngularity_R{jetR}").Fill(value)

                # Jet mass
                # NOTE: Charged jet pt
                if 60 < jet_pt < 80:
                    # Could also use modp2(), but it would make it asymmetric compared to the squared E,
                    # which will for sure make me do a double take at some point. Better to avoid it.
                    jet_mass = np.sqrt(jet.E() ** 2 - jet.pt() ** 2 - jet.pz() ** 2)

                    # Add holes together as four vectors, and then calculate their mass, removing it from the jet mass.
                    if holes_in_jet:
                        # Copy to avoid modifying the original hole.
                        hole_four_vector = fj.PseudoJet(holes_in_jet[0])
                        for hadron in holes_in_jet[1:]:
                            hole_four_vector += hadron

                        # Remove mass from holes
                        jet_mass -= np.sqrt(hole_four_vector.E() ** 2 - hole_four_vector.pt() ** 2 - hole_four_vector.pz() ** 2)

                    getattr(self, f"hJetMass_R{jetR}").Fill(jet_mass)

                # Substructure
                # NOTE: Charged jet pt
                if 60 < jet_pt < 80:
                    reclustering_algorithm = fj.cambridge_algorithm
                    gshop = fjcontrib.GroomerShop(jet, jetR, reclustering_algorithm)
                    jet_groomed_lund = gshop.soft_drop(self.soft_drop_beta, self.soft_drop_zcut, jetR)
                    # TODO: Careful with untagged bin...
                    if jet_groomed_lund:
                        theta_g = jet_groomed_lund.Delta() / jetR
                        zg = jet_groomed_lund.z()
                    else:
                        # Following my own convention. Revise as necessary.
                        theta_g = -0.05
                        zg = -0.05
                    getattr(self, f"hSoftDrop_zg_R{jetR}").Fill(zg)
                    getattr(self, f"hSoftDrop_theta_g_R{jetR}").Fill(theta_g)

                # h-jet
                # TODO: Caution - this for sure needs a second set of eyes.
                # NOTE: We may need a lower jet pt range to determine cref, but I didn't look into this.
                for hadron in itertools.chain(fj_hadrons_positive, fj_hadrons_negative):
                    found_low = False
                    found_high = False
                    # TODO: Doesn't account for detector effects on hadron.
                    if self.hjet_low_trigger_range[0] < hadron.pt() < self.hjet_low_trigger_range[1]:
                        found_low = True
                    if self.hjet_high_trigger_range[0] < hadron.pt() < self.hjet_high_trigger_range[1]:
                        found_high = True

                    # Check if jet is within the right range.
                    found_recoil_jet = False
                    if found_low or found_high:
                        if np.pi - jet.delta_R(hadron) < 0.6:
                            found_recoil_jet = True

                    # Record the jet pt if it was found
                    if found_recoil_jet:
                        # Jet yield
                        hist_label = "high" if found_high else "low"
                        getattr(self, f"hRecoilJetYield_{hist_label}Trigger").Fill(jet_pt)
                        # Delta Phi
                        # NOTE: Charged jet pt
                        if 40 < jet_pt < 60:
                            getattr(self, f"hRecoilJetDeltaPhi_{hist_label}Trigger").Fill(hadron.delta_phi_to(jet))

            # Recoil histogram
            getattr(self, 'hJetPt_recoils_R{}'.format(jetR)).Fill(jet_pt, negative_pt)


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

    analysis = AnalyzeJetscapeEvents_PHYS(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
