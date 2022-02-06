#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE parquet output file,
  and write out a new parquet file containing calculated observables

  For AA, must perform hole subtraction:
    For hadron observables:
      We save spectra for positive/negative particles separately, then subtract at histogram-level in plotting script
    For jets we find three different collections of jets:
      (1) Using shower+recoil particles, with constituent subtraction
           - No further hole subtraction necessary
      (2) Using shower+recoil particles, with standard recombiner
          In this case, observable-specific hole subtraction necessary
          We consider three different classes of jet observables:
           (i) Jet pt-like observables -- subtract holes within R
           (ii) Additive substructure -- subtract holes within R
           (iii) Non-additive substructure -- correct the jet pt only
          We also save unsubtracted histograms for comparison (although for substructure we still correct pt)
      (3) Using shower+recoil+hole particles, with negative recombiner
          In this case, observable-specific hole subtraction necessary
          We consider three different classes of jet observables:
           (i) Jet pt-like observables -- no further hole subtraction
           (ii) Additive substructure -- subtract holes within R
           (iii) Non-additive substructure -- no further hole subtraction

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
from collections import defaultdict
from pathlib import Path
from collections import defaultdict

# Fastjet via python (from external library heppy)
import fjcontrib
import fjext
import fastjet as fj

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
        self.output_file = f'observables'
        # Update the output_file to contain the labeling in the final_state_hadrons file.
        # We use this naming convention as the flag for whether we should attempt to rename it.
        if "final_state_hadrons" in self.input_file_hadrons:
            _input_filename = Path(self.input_file_hadrons).name
            # The filename will be something like "observables_0000_00.parquet", assuming
            # that the original name was "observables"
            self.output_file = _input_filename.replace("final_state_hadrons", self.output_file)
            #print(f'Updated output_file name to "{self.output_file}" in order to add identifying indices.')

        # Load observable blocks
        self.hadron_observables = config['hadron']
        self.hadron_correlation_observables = config['hadron_correlations']
        self.inclusive_chjet_observables = config['inclusive_chjet']
        self.inclusive_jet_observables = {}
        self.semi_inclusive_chjet_observables = {}
        self.dijet_observables = {}
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

        # General grooming parameters'
        self.grooming_settings = {}
        if 'SoftDrop' in config:
            self.grooming_settings = config['SoftDrop']

        # If AA, set different options for hole subtraction treatment
        if self.is_AA:
            self.jet_collection_labels = config['jet_collection_labels']
        else:
            self.jet_collection_labels = ['']

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Initialize a dictionary that will store a list of calculated values for each output observable
        self.observable_dict_event = defaultdict(list)

        # Create list of fastjet::PseudoJets (separately for jet shower particles and holes)
        fj_hadrons_positive, pid_hadrons_positive = self.fill_fastjet_constituents(event, select_status='+')
        fj_hadrons_negative, pid_hadrons_negative = self.fill_fastjet_constituents(event, select_status='-')

        # Create list of charged particles
        fj_hadrons_positive_charged, pid_hadrons_positive_charged = self.fill_fastjet_constituents(event, select_status='+',
                                                                     select_charged=True)
        fj_hadrons_negative_charged, pid_hadrons_negative_charged = self.fill_fastjet_constituents(event, select_status='-',
                                                                     select_charged=True)

        # Fill hadron observables for jet shower particles
        self.fill_hadron_observables(fj_hadrons_positive, pid_hadrons_positive, status='+')
        if self.is_AA:
            self.fill_hadron_observables(fj_hadrons_negative, pid_hadrons_negative, status='-')

        # Fill hadron correlation observables
        self.fill_hadron_correlation_observables(fj_hadrons_positive, pid_hadrons_positive)

        # Fill jet observables
        for jet_collection_label in self.jet_collection_labels:

            # If constituent subtraction, subtract the event (with rho determined from holes) -- we can then neglect the holes
            if jet_collection_label == '_constituent_subtraction':
                self.bge_rho.set_particles(fj_hadrons_negative)
                hadrons_positive = self.constituent_subtractor.subtract_event(fj_hadrons_positive)
                hadrons_negative = None

                self.bge_rho.set_particles(fj_hadrons_negative_charged)
                hadrons_positive_charged = self.constituent_subtractor.subtract_event(fj_hadrons_positive_charged)
                hadrons_negative_charged = None

            # For shower_recoil and negative_recombiner cases, keep both positive and negative hadrons
            else:
                hadrons_positive = fj_hadrons_positive
                hadrons_negative = fj_hadrons_negative
                hadrons_positive_charged = fj_hadrons_positive_charged
                hadrons_negative_charged = fj_hadrons_negative_charged

            # Find jets and fill observables
            self.fill_jet_observables(hadrons_positive, hadrons_negative,
                                      hadrons_positive_charged, hadrons_negative_charged,
                                      pid_hadrons_positive, pid_hadrons_negative,
                                      pid_hadrons_positive_charged, pid_hadrons_negative_charged,
                                      jet_collection_label=jet_collection_label)

    # ---------------------------------------------------------------
    # Fill hadron observables
    # (assuming weak strange decays are off, but charm decays are on)
    # ---------------------------------------------------------------
    def fill_hadron_observables(self, fj_particles, pid_hadrons, status='+'):

        # Note that for identified particles, we store holes of the identified species
        suffix = ''
        if status == '-':
            suffix = '_holes'

        # Loop through hadrons
        for particle in fj_particles:

            # Fill some basic hadron info
            pid = pid_hadrons[np.abs(particle.user_index())-1]
            pt = particle.pt()
            eta = particle.eta()

            if self.sqrts in [2760, 5020]:

                # ALICE
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if self.centrality_accepted(self.hadron_observables['pt_ch_alice']['centrality']):
                    pt_min = self.hadron_observables['pt_ch_alice']['pt'][0]
                    pt_max = self.hadron_observables['pt_ch_alice']['pt'][1]
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_ch_alice']['eta_cut']:
                            if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                self.observable_dict_event[f'hadron_pt_ch_alice{suffix}'].append(pt)

                # Charged pion
                if self.centrality_accepted(self.hadron_observables['pt_pi_alice']['centrality']):
                    pt_min = self.hadron_observables['pt_pi_alice']['pt'][0]
                    pt_max = self.hadron_observables['pt_pi_alice']['pt'][1]
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_pi_alice']['eta_cut']:
                            if abs(pid) == 211:
                                self.observable_dict_event[f'hadron_pt_pi_alice{suffix}'].append(pt)

                # Neutral pions
                if self.sqrts in [2760]:
                    if self.centrality_accepted(self.hadron_observables['pt_pi0_alice']['centrality']):
                        pt_min = self.hadron_observables['pt_pi0_alice']['pt'][0]
                        pt_max = self.hadron_observables['pt_pi0_alice']['pt'][1]
                        if pt > pt_min and pt < pt_max:
                            if abs(eta) < self.hadron_observables['pt_pi0_alice']['eta_cut']:
                                if abs(pid) == 111:
                                    self.observable_dict_event[f'hadron_pt_pi0_alice{suffix}'].append(pt)

                # ATLAS
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if self.sqrts in [2760]:
                    if self.centrality_accepted(self.hadron_observables['pt_ch_atlas']['centrality']):
                        pt_min = self.hadron_observables['pt_ch_atlas']['pt'][0]
                        pt_max = self.hadron_observables['pt_ch_atlas']['pt'][1]
                        if pt > pt_min and pt < pt_max:
                            if abs(eta) < self.hadron_observables['pt_ch_atlas']['eta_cut']:
                                if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                    self.observable_dict_event[f'hadron_pt_ch_atlas{suffix}'].append(pt)

                # CMS
                # Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if self.centrality_accepted(self.hadron_observables['pt_ch_cms']['centrality']):
                    pt_min = self.hadron_observables['pt_ch_cms']['pt'][0]
                    pt_max = self.hadron_observables['pt_ch_cms']['pt'][1]
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_ch_cms']['eta_cut']:
                            if abs(pid) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                self.observable_dict_event[f'hadron_pt_ch_cms{suffix}'].append(pt)

            elif self.sqrts in [200]:

                # PHENIX
                # Neutral pions
                if self.centrality_accepted(self.hadron_observables['pt_pi0_phenix']['centrality']):
                    pt_min = self.hadron_observables['pt_pi0_phenix']['pt'][0]
                    pt_max = 100. # Open upper bound
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_pi0_phenix']['eta_cut']:
                            if abs(pid) == 111:
                                self.observable_dict_event[f'hadron_pt_pi0_phenix{suffix}'].append(pt)

                # STAR
                # Charged hadrons (pi+, K+, p+)
                if self.centrality_accepted(self.hadron_observables['pt_ch_star']['centrality']):
                    pt_min = self.hadron_observables['pt_ch_star']['pt'][0]
                    pt_max = 100. # Open upper bound
                    if pt > pt_min and pt < pt_max:
                        if abs(eta) < self.hadron_observables['pt_ch_star']['eta_cut']:
                            if abs(pid) in [211, 321, 2212]:
                                self.observable_dict_event[f'hadron_pt_ch_star{suffix}'].append(pt)

    # ---------------------------------------------------------------
    # Fill hadron correlation observables
    # ---------------------------------------------------------------
    def fill_hadron_correlation_observables(self, fj_particles, pid_hadrons) -> None:
        # NOTE: The loop order here is different than other functions because without some optimization,
        #       it's very easy to have an O(n^2) loop looking for trigger and associated particles.
        #       We keep track of the particles which pass our conditions, so the double loop ends up
        #       running for small n rather than the full event multiplicity
        if self.sqrts in [200]:
            if self.centrality_accepted(self.hadron_correlation_observables['dihadron_star']['centrality']):  # type: ignore
                # Keep track of the trigger particles in pt ranges
                trigger_particles = defaultdict(list)
                # Keep track of the associated particles in pt ranges
                associated_particles = defaultdict(list)
                pt_trigger_ranges = self.hadron_correlation_observables["dihadron_star"]["pt_trig"]
                pt_associated_ranges = self.hadron_correlation_observables["dihadron_star"]["pt_assoc"]
                for i, particle in enumerate(fj_particles):
                    # eta cut
                    if abs(particle.eta()) < self.hadron_correlation_observables['dihadron_star']['eta_cut']:
                        # Charged hadrons (pi+, K+, p+)
                        if pid_hadrons[np.abs(particle.user_index())-1] in [211, 321, 2212]:
                            pt = particle.pt()
                            for pt_trig_range in pt_trigger_ranges:
                                pt_trig_min, pt_trig_max = pt_trig_range
                                if pt_trig_min <= pt < pt_trig_max:
                                    # Found trigger - save it
                                    trigger_particles[(pt_trig_min, pt_trig_max)].append(particle)

                            for pt_assoc_range in pt_associated_ranges:
                                pt_assoc_min, pt_assoc_max = pt_assoc_range
                                # If the upper range has -1, it's unbounded, so we make it large enough not to matter
                                pt_assoc_max = 1000 if pt_assoc_max == -1 else pt_assoc_max
                                if pt_assoc_min <= pt < pt_assoc_max:
                                    associated_particles[(pt_assoc_min, pt_assoc_max)].append(particle)

                # Now, create the correlations over our reduced set of particles
                for (pt_trig_min, pt_trig_max), trig_particles in trigger_particles.items():
                    for trigger_particle in trig_particles:
                        # Store the trigger pt to count the number of triggers
                        # In principle, the dphi list could have been enough to get the number of triggers,
                        # but it's not so easy to integrate with the existing histogram code.
                        # So we keep separate track of the triggers, same as is done for D(z)
                        self.observable_dict_event[f'hadron_correlations_dihadron_star_Ntrig'].append(trigger_particle.pt())

                        for (pt_assoc_min, pt_assoc_max), assoc_particles in associated_particles.items():
                            # First, just calculate the values
                            dphi_values = []
                            for associated_particle in assoc_particles:
                                # Trigger particle pt must be larger than the associated particle.
                                # If it's smaller, skip it. We'll have picked up the higher pt particle as a trigger, so
                                # we'll account for it when we loop over that trigger.
                                # NOTE: We also want to ensure that the trigger and associated aren't the same particle.
                                #       However, requiring trigger pt > assoc pt also implicitly requires that the
                                #       two particles can't be the same.
                                if trigger_particle.pt() > associated_particle.pt():
                                    # stores phi_trig - phi_assoc
                                    dphi_values.append(associated_particle.delta_phi_to(trigger_particle))

                            # If nothing to record, then skip over this trigger
                            # NOTE: This is actually quite important for the output because calling extend on a defaultdict
                            #       with an empty list causes the None (what we want for an empty list for our output) to
                            #       be stored as an empty list (which we don't want). Unfortunately, this is subtle to see
                            #       when debugging, so have some care and keep an eye out for this with this observable.
                            if len(dphi_values) == 0:
                                continue

                            # Label with both pt ranges
                            label = f"pt_trig_{pt_trig_min:g}_{pt_trig_max:g}_pt_assoc_{pt_assoc_min:g}_{pt_assoc_max:g}"
                            # Store a list of dphi of associated particles
                            # NOTE: We actually store the triggers separately, but in principle we could extract it directly from here if
                            #       we stored the associated particles per trigger (ie. used append rather than extend). However,
                            #       it's more difficult to integrate with the existing infrastructure, but so we take the simpler route
                            #       and use a flat list
                            # NOTE: Here we standardize the values to match with the measured correlation range
                            self.observable_dict_event[f'hadron_correlations_dihadron_star_{label}'].extend([
                                analyze_events_base_STAT.dphi_in_range_for_hadron_correlations(phi) for phi in dphi_values
                            ])

    # ---------------------------------------------------------------
    # Fill jet observables
    # For AA, we find three different collections of jets:
    #
    #   (1) Using shower+recoil particles, with constituent subtraction
    #        - No further hole subtraction necessary
    #
    #   (2) Using shower+recoil particles, using standard recombiner
    #       In this case, observable-specific hole subtraction necessary
    #       We consider three different classes of jet observables:
    #        (i) Jet pt-like observables -- subtract holes within R
    #        (ii) Additive substructure -- subtract holes within R
    #        (iii) Non-additive substructure -- correct the jet pt only
    #       We also save unsubtracted histograms for comparison.
    #
    #   (3) Using shower+recoil+hole particles, using negative recombiner
    #       In this case, observable-specific hole subtraction necessary
    #       We consider three different classes of jet observables:
    #        (i) Jet pt-like observables -- no further hole subtraction
    #        (ii) Additive substructure -- subtract holes within R
    #        (iii) Non-additive substructure -- we do no further hole subtraction
    # ---------------------------------------------------------------
    def fill_jet_observables(self, hadrons_positive, hadrons_negative,
                             hadrons_positive_charged, hadrons_negative_charged,
                             pid_hadrons_positive, pid_hadrons_negative,
                             pid_hadrons_positive_charged, pid_hadrons_negative_charged,
                             jet_collection_label=''):

        # Set the appropriate lists of hadrons to input to the jet finding
        if jet_collection_label in ['', '_shower_recoil', '_constituent_subtraction']:
            hadrons_for_jet_finding = hadrons_positive
            hadrons_for_jet_finding_charged = hadrons_positive_charged
        elif jet_collection_label in ['_negative_recombiner']:
            hadrons_for_jet_finding = list(hadrons_positive) + list(hadrons_negative)
            hadrons_for_jet_finding_charged = list(hadrons_positive_charged) + list(hadrons_negative_charged)

        # Loop through specified jet R
        for jetR in self.jet_R:

            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            if jet_collection_label in ['_negative_recombiner']:
                recombiner = fjext.NegativeEnergyRecombiner()
                jet_def.set_recombiner(recombiner)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(self.max_jet_y)

            # Full jets
            self.find_jets_and_fill(hadrons_for_jet_finding, hadrons_negative,
                                    pid_hadrons_positive, pid_hadrons_negative,
                                    jet_def, jet_selector, jetR, jet_collection_label, full_jet=True)

            # Charged jets
            self.find_jets_and_fill(hadrons_for_jet_finding_charged, hadrons_negative_charged,
                                    pid_hadrons_positive_charged, pid_hadrons_negative_charged,
                                    jet_def, jet_selector, jetR, jet_collection_label, full_jet=False)

    # ---------------------------------------------------------------
    # Find jets and fill histograms -- either full or charged
    # ---------------------------------------------------------------
    def find_jets_and_fill(self, hadrons_for_jet_finding, hadrons_negative,
                           pid_hadrons_positive, pid_hadrons_negative,
                           jet_def, jet_selector, jetR, jet_collection_label, full_jet=True):

        # Fill inclusive jets
        cs = fj.ClusterSequence(hadrons_for_jet_finding, jet_def)
        jets = fj.sorted_by_pt(cs.inclusive_jets())
        jets_selected = jet_selector(jets)

        [self.analyze_inclusive_jet(jet, hadrons_for_jet_finding, hadrons_negative,
                                    pid_hadrons_positive, pid_hadrons_negative,
                                    jetR, full_jet=full_jet,
                                    jet_collection_label=jet_collection_label) for jet in jets_selected]

        # Fill dijet observables -- full jets only
        if full_jet:
            if self.dijet_observables:
                self.fill_dijet_observables(jets_selected, hadrons_negative, jetR,
                                            jet_collection_label=jet_collection_label)

        # Fill semi-inclusive jet correlations -- charged jets only
        if not full_jet:
            if self.semi_inclusive_chjet_observables:
                if self.sqrts == 2760:
                    jetR_list = self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']+self.semi_inclusive_chjet_observables['nsubjettiness_alice']['jet_R']
                elif self.sqrts == 200:
                    jetR_list = self.semi_inclusive_chjet_observables['IAA_star']['jet_R']
                if jetR in jetR_list:
                    self.fill_semi_inclusive_chjet_observables(jets_selected, hadrons_for_jet_finding, hadrons_negative,
                                                                jetR, jet_collection_label=jet_collection_label)

    # ---------------------------------------------------------------
    # Fill inclusive jet observables
    # ---------------------------------------------------------------
    def analyze_inclusive_jet(self, jet, hadrons_for_jet_finding, hadrons_negative,
                              pid_hadrons_positive, pid_hadrons_negative,
                              jetR, full_jet=True, jet_collection_label=''):

        # Get the list of holes inside the jet, if applicable
        #   For the shower+recoil case, we need to subtract the hole pt
        #   For the negative recombiner case, we do not need to adjust the pt, but we want to keep track of the holes
        holes_in_jet = []
        if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
            for hadron in hadrons_negative:
                if jet.delta_R(hadron) < jetR:
                    holes_in_jet.append(hadron)

        # Correct the pt of the jet, if applicable
        # For pp or negative recombiner or constituent subtraction case, we do not need to adjust the pt
        # For the shower+recoil case, we need to subtract the hole pt
        if jet_collection_label in ['', '_negative_recombiner', '_constituent_subtraction']:
            jet_pt = jet_pt_uncorrected = jet.pt()
        elif jet_collection_label in ['_shower_recoil']:
            negative_pt = 0.
            for hadron in holes_in_jet:
                negative_pt += hadron.pt()
            jet_pt_uncorrected = jet.pt()               # uncorrected pt: shower+recoil
            jet_pt = jet_pt_uncorrected - negative_pt   # corrected pt: shower+recoil-holes

        # Fill observables
        if full_jet:

            # Ungroomed
            self.fill_full_jet_ungroomed_observables(jet, hadrons_for_jet_finding, holes_in_jet,
                                                     pid_hadrons_positive, pid_hadrons_negative,
                                                     jet_pt, jet_pt_uncorrected, jetR, jet_collection_label=jet_collection_label)

            # Groomed
            if self.grooming_settings:
                for grooming_setting in self.grooming_settings:
                    self.fill_full_jet_groomed_observables(grooming_setting, jet, jet_pt, jetR, jet_collection_label=jet_collection_label)

        else:

            # Ungroomed
            self.fill_charged_jet_ungroomed_observables(jet, holes_in_jet, pid_hadrons_positive, jet_pt, jet_pt_uncorrected,
                                                        jetR, jet_collection_label=jet_collection_label)

            # Groomed
            if self.grooming_settings:
                for grooming_setting in self.grooming_settings:
                    self.fill_charged_jet_groomed_observables(grooming_setting, jet, jet_pt, jetR,
                                                              jet_collection_label=jet_collection_label)

    # ---------------------------------------------------------------
    # Fill inclusive full jet observables
    # ---------------------------------------------------------------
    def fill_full_jet_ungroomed_observables(self, jet, hadrons_for_jet_finding, holes_in_jet,
                                            pid_hadrons_positive, pid_hadrons_negative,
                                            jet_pt, jet_pt_uncorrected, jetR, jet_collection_label=''):

        if self.sqrts in [2760, 5020]:

            # ALICE RAA
            #   Hole treatment:
            #    - For RAA, all jet collections can be filled from the corrected jet pt
            #    - In the shower_recoil case, we also fill the unsubtracted jet pt
            if self.centrality_accepted(self.inclusive_jet_observables['pt_alice']['centrality']):
                pt_min, pt_max = self.inclusive_jet_observables['pt_alice']['pt']
                if jetR in self.inclusive_jet_observables['pt_alice']['jet_R']:
                    if abs(jet.eta()) < (self.inclusive_jet_observables['pt_alice']['eta_cut_R'] - jetR):
                        if pt_min < jet_pt < pt_max:

                            # Check leading track requirement
                            if jetR == 0.2:
                                min_leading_track_pt = 5.
                            else:
                                min_leading_track_pt = 7.

                            accept_jet = False
                            acceptable_hadrons = [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]
                            for constituent in jet.constituents():
                                if constituent.pt() > min_leading_track_pt:
                                    # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                    if abs(pid_hadrons_positive[np.abs(constituent.user_index())-1]) in acceptable_hadrons:
                                        accept_jet = True
                            if accept_jet:
                                self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}{jet_collection_label}'].append(jet_pt)
                                if jet_collection_label in ['_shower_recoil']:
                                    self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

            # ATLAS RAA
            if self.centrality_accepted(self.inclusive_jet_observables['pt_atlas']['centrality']):
                pt_min = self.inclusive_jet_observables['pt_atlas']['pt'][0]
                pt_max = self.inclusive_jet_observables['pt_atlas']['pt'][1]
                if jetR in self.inclusive_jet_observables['pt_atlas']['jet_R']:
                    if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
                        if pt_min < jet_pt < pt_max:
                            self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}{jet_collection_label}'].append(jet_pt)
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

            # ATLAS RAA -- rapidity-dependence
            if self.centrality_accepted(self.inclusive_jet_observables['pt_y_atlas']['centrality']):
                pt_min = self.inclusive_jet_observables['pt_y_atlas']['pt'][0]
                pt_max = self.inclusive_jet_observables['pt_y_atlas']['pt'][-1]
                if jetR in self.inclusive_jet_observables['pt_y_atlas']['jet_R']:
                    y_abs = abs(jet.rap())
                    if y_abs < self.inclusive_jet_observables['pt_y_atlas']['y_cut']:
                        if pt_min < jet_pt < pt_max:
                            self.observable_dict_event[f'inclusive_jet_pt_y_atlas_R{jetR}{jet_collection_label}'].append([jet_pt, y_abs])
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_jet_pt_y_atlas_R{jetR}{jet_collection_label}_unsubtracted'].append([jet_pt_uncorrected, y_abs])

            # CMS RAA
            if self.centrality_accepted(self.inclusive_jet_observables['pt_cms']['centrality']):
                pt_min = self.inclusive_jet_observables['pt_cms']['pt'][0]
                pt_max = self.inclusive_jet_observables['pt_cms']['pt'][1]
                if jetR in self.inclusive_jet_observables['pt_cms']['jet_R']:
                    if abs(jet.eta()) < self.inclusive_jet_observables['pt_cms']['eta_cut']:
                        if pt_min < jet_pt < pt_max:
                            self.observable_dict_event[f'inclusive_jet_pt_cms_R{jetR}{jet_collection_label}'].append(jet_pt)
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_jet_pt_cms_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

            # ATLAS D(z)
            #   Hole treatment:
            #    - For show_recoil case, we separately fill using hadrons_for_jet_finding (which are positive only) and holes_in_jet
            #    - For negative_recombiner case, we separately fill the positive-status and negative-status hadrons_for_jet_finding
            #    - For constituent_subtraction, we will using hadrons_for_jet_finding (which are positive only)
            #   Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            acceptable_hadrons = [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]
            if self.centrality_accepted(self.inclusive_jet_observables['Dz_atlas']['centrality']):
                pt_min = self.inclusive_jet_observables['Dz_atlas']['pt'][0]
                pt_max = self.inclusive_jet_observables['Dz_atlas']['pt'][-1]
                if jetR in self.inclusive_jet_observables['Dz_atlas']['jet_R']:
                    if abs(jet.rap()) < self.inclusive_jet_observables['Dz_atlas']['y_cut']:
                        if pt_min < jet_pt < pt_max:
                            self.observable_dict_event[f'inclusive_jet_Dz_atlas_R{jetR}{jet_collection_label}_Njets'].append(jet_pt)
                            for hadron in hadrons_for_jet_finding:
                                if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() < 0 :
                                    continue
                                pid = pid_hadrons_positive[np.abs(hadron.user_index())-1]
                                if abs(pid) in acceptable_hadrons:
                                    if jet.delta_R(hadron) < jetR:
                                        z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                        self.observable_dict_event[f'inclusive_jet_Dz_atlas_R{jetR}{jet_collection_label}'].append([jet_pt, z])
                                        self.observable_dict_event[f'inclusive_jet_Dpt_atlas_R{jetR}{jet_collection_label}'].append([jet_pt, hadron.pt()])
                            if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                for hadron in holes_in_jet:
                                    if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                        continue
                                    pid = pid_hadrons_negative[np.abs(hadron.user_index())-1]
                                    if abs(pid) in acceptable_hadrons:
                                        if jet.delta_R(hadron) < jetR:
                                            z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                            self.observable_dict_event[f'inclusive_jet_Dz_atlas_R{jetR}_holes{jet_collection_label}'].append([jet_pt, z])
                                            self.observable_dict_event[f'inclusive_jet_Dpt_atlas_R{jetR}_holes{jet_collection_label}'].append([jet_pt, hadron.pt()])

            # CMS D(z)
            #   Hole treatment:
            #    - For show_recoil case, we separately fill using hadrons_for_jet_finding (which are positive only) and holes_in_jet
            #    - For negative_recombiner case, we separately fill the positive-status and negative-status hadrons_for_jet_finding
            #    - For constituent_subtraction, we will using hadrons_for_jet_finding (which are positive only)
            #   Charged hadrons (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            acceptable_hadrons = [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]
            if self.sqrts == 2760:
                if self.centrality_accepted(self.inclusive_jet_observables['Dz_cms']['centrality']):
                    pt_min = self.inclusive_jet_observables['Dz_cms']['pt'][0]
                    pt_max = self.inclusive_jet_observables['Dz_cms']['pt'][-1]
                    eta_range = self.inclusive_jet_observables['Dz_cms']['eta_cut']
                    track_pt_min = self.inclusive_jet_observables['Dz_cms']['track_pt_min']
                    if jetR in self.inclusive_jet_observables['Dz_cms']['jet_R']:
                        if eta_range[0] < abs(jet.eta()) < eta_range[1]:
                            if pt_min < jet_pt < pt_max:
                                self.observable_dict_event[f'inclusive_jet_Dz_cms_R{jetR}{jet_collection_label}_Njets'].append(jet_pt)
                                for hadron in hadrons_for_jet_finding:
                                    if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() < 0 :
                                        continue
                                    if hadron.pt() > track_pt_min:
                                        pid = pid_hadrons_positive[np.abs(hadron.user_index())-1]
                                        if abs(pid) in acceptable_hadrons:
                                            if jet.delta_R(hadron) < jetR:
                                                z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                                xi = np.log(1/z)
                                                self.observable_dict_event[f'inclusive_jet_Dz_cms_R{jetR}{jet_collection_label}'].append([jet_pt, xi])
                                                self.observable_dict_event[f'inclusive_jet_Dpt_cms_R{jetR}{jet_collection_label}'].append([jet_pt, hadron.pt()])
                                if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                    for hadron in holes_in_jet:
                                        if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                            continue
                                        if hadron.pt() > track_pt_min:
                                            pid = pid_hadrons_negative[np.abs(hadron.user_index())-1]
                                            if abs(pid) in acceptable_hadrons:
                                                if jet.delta_R(hadron) < jetR:
                                                    z = hadron.pt() * np.cos(jet.delta_R(hadron)) / jet_pt
                                                    xi = np.log(1/z)
                                                    self.observable_dict_event[f'inclusive_jet_Dz_cms_R{jetR}_holes{jet_collection_label}'].append([jet_pt, xi])
                                                    self.observable_dict_event[f'inclusive_jet_Dpt_cms_R{jetR}_holes{jet_collection_label}'].append([jet_pt, hadron.pt()])

            # CMS jet charge
            #   Hole treatment:
            #    - For show_recoil case, we subtract the contribution of holes within R (and also store the unsubtracted charge)
            #    - For negative_recombiner case, we subtract the contribution of holes within R
            #    - For constituent_subtraction, no subtraction is needed
            # Charged particles (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            acceptable_hadrons = [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]
            if self.sqrts == 5020:
                if self.centrality_accepted(self.inclusive_jet_observables['charge_cms']['centrality']):
                    pt_min = self.inclusive_jet_observables['charge_cms']['pt_min']
                    if jetR in self.inclusive_jet_observables['charge_cms']['jet_R']:
                        if abs(jet.eta()) < self.inclusive_jet_observables['charge_cms']['eta_cut']:
                            if jet_pt > pt_min:
                                for kappa in self.inclusive_jet_observables['charge_cms']['kappa']:
                                    sum = 0
                                    for hadron in hadrons_for_jet_finding:
                                        if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() < 0 :
                                            continue
                                        if hadron.pt() > self.inclusive_jet_observables['charge_cms']['track_pt_min']:
                                            pid = pid_hadrons_positive[np.abs(hadron.user_index())-1]
                                            if abs(pid) in acceptable_hadrons:
                                                if jet.delta_R(hadron) < jetR:
                                                    sum += self.charge(pid) * np.power(hadron.pt(), kappa)
                                    if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                        sum_holes = 0
                                        for hadron in holes_in_jet:
                                            if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                                continue
                                            if hadron.pt() > self.inclusive_jet_observables['charge_cms']['track_pt_min']:
                                                pid = pid_hadrons_negative[np.abs(hadron.user_index())-1]
                                                if abs(pid) in acceptable_hadrons:
                                                    if jet.delta_R(hadron) < jetR:
                                                        sum_holes += self.charge(pid) * np.power(hadron.pt(), kappa)
                                        charge = (sum - sum_holes) / np.power(jet_pt, kappa)
                                        if jet_collection_label in ['_shower_recoil']:
                                            charge_unsubtracted = sum / np.power(jet_pt, kappa)
                                            self.observable_dict_event[f'inclusive_jet_charge_cms_R{jetR}_k{kappa}{jet_collection_label}_unsubtracted'].append(charge_unsubtracted)
                                    else:
                                        charge = sum / np.power(jet_pt, kappa)
                                    self.observable_dict_event[f'inclusive_jet_charge_cms_R{jetR}_k{kappa}{jet_collection_label}'].append(charge)

    # ---------------------------------------------------------------
    # Fill inclusive full jet observables
    # ---------------------------------------------------------------
    def fill_full_jet_groomed_observables(self, grooming_setting, jet, jet_pt, jetR, jet_collection_label=''):

        # Construct groomed jet

        # For negative_recombiner case, we set the negative recombiner also for the C/A reclustering
        jet_def = fj.JetDefinition(fj.cambridge_algorithm, jetR)
        if jet_collection_label in ['_negative_recombiner']:
            recombiner = fjext.NegativeEnergyRecombiner()
            jet_def.set_recombiner(recombiner)
        gshop = fjcontrib.GroomerShop(jet, jet_def)

        zcut = grooming_setting['zcut']
        beta = grooming_setting['beta']
        jet_groomed_lund = gshop.soft_drop(beta, zcut, jetR)
        if not jet_groomed_lund:
            return

        if self.sqrts == 5020:

            # CMS m_g
            #   Hole treatment:
            #    - For show_recoil case, correct the pt only
            #    - For negative_recombiner case, no subtraction is needed
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_jet_observables['mg_cms']['centrality']):
                if grooming_setting in self.inclusive_jet_observables['mg_cms']['SoftDrop']:
                    pt_min = self.inclusive_jet_observables['mg_cms']['pt'][0]
                    pt_max = self.inclusive_jet_observables['mg_cms']['pt'][-1]
                    if jetR in self.inclusive_jet_observables['mg_cms']['jet_R']:
                        if abs(jet.eta()) < (self.inclusive_jet_observables['mg_cms']['eta_cut']):
                            if pt_min < jet_pt < pt_max:
                                if jet_groomed_lund.Delta() > self.inclusive_jet_observables['mg_cms']['dR']:
                                    mg = jet_groomed_lund.pair().m() / jet_pt    # Note: untagged jets will return negative value
                                    self.observable_dict_event[f'inclusive_jet_mg_cms_R{jetR}_zcut{zcut}_beta{beta}{jet_collection_label}'].append([jet_pt, mg])

            # CMS z_g
            #   Hole treatment:
            #    - For show_recoil case, correct the pt only
            #    - For negative_recombiner case, no subtraction is needed
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_jet_observables['zg_cms']['centrality']):
                if grooming_setting in self.inclusive_jet_observables['zg_cms']['SoftDrop']:
                    pt_min = self.inclusive_jet_observables['zg_cms']['pt'][0]
                    pt_max = self.inclusive_jet_observables['zg_cms']['pt'][-1]
                    if jetR in self.inclusive_jet_observables['zg_cms']['jet_R']:
                        if abs(jet.eta()) < (self.inclusive_jet_observables['zg_cms']['eta_cut']):
                            if pt_min < jet_pt < pt_max:
                                if jet_groomed_lund.Delta() > self.inclusive_jet_observables['zg_cms']['dR']:
                                    zg = jet_groomed_lund.z()    # Note: untagged jets will return negative value
                                    self.observable_dict_event[f'inclusive_jet_zg_cms_R{jetR}_zcut{zcut}_beta{beta}{jet_collection_label}'].append([jet_pt, zg])

    # ---------------------------------------------------------------
    # Fill inclusive charged jet observables
    # ---------------------------------------------------------------
    def fill_charged_jet_ungroomed_observables(self, jet, holes_in_jet, pid_hadrons_positive,
                                               jet_pt, jet_pt_uncorrected, jetR, jet_collection_label=''):

        if self.sqrts == 5020:

            # ALICE subjet z_R
            #   Hole treatment:
            #    - For show_recoil case, subtract holes within r (for subjets) and R (for jets)
            #    - For negative_recombiner case, subtract holes within r (for subjets) only
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['zr_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['zr_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['zr_alice']['pt'][-1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['zr_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['zr_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:
                            for r in self.inclusive_chjet_observables['zr_alice']['r']:

                                cs_subjet = fj.ClusterSequence(jet.constituents(), fj.JetDefinition(fj.antikt_algorithm, r))
                                subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
                                _, leading_subjet_pt, _ = self.leading_jet(subjets, holes_in_jet, r)
                                z_leading = leading_subjet_pt / jet_pt
                                if np.isclose(z_leading, 1.): # If z=1, it will be default be placed in overflow bin -- prevent this
                                    z_leading = 0.999

                                self.observable_dict_event[f'inclusive_chjet_zr_alice_R{jetR}_r{r}{jet_collection_label}'].append([z_leading])

            # ALICE jet axis Standard-WTA
            #   Hole treatment:
            #    - For show_recoil case, correct the pt only
            #    - For negative_recombiner case, no subtraction is needed, although we recluster using the negative recombiner again
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['axis_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['axis_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['axis_alice']['pt'][-1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['axis_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['axis_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:

                            jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
                            if jet_collection_label in ['_negative_recombiner']:
                                recombiner = fjext.NegativeEnergyRecombiner()
                                jet_def_wta.set_recombiner(recombiner)
                            jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
                            reclusterer_wta = fjcontrib.Recluster(jet_def_wta)
                            jet_wta = reclusterer_wta.result(jet)
                            deltaR = jet.delta_R(jet_wta)

                            self.observable_dict_event[f'inclusive_chjet_axis_alice_R{jetR}{jet_collection_label}'].append([jet_pt, deltaR])

            # ALICE ungroomed angularity
            #   Hole treatment:
            #    - For show_recoil case, subtract the hole contribution within R to the angularity (also store unsubtracted case)
            #    - For negative_recombiner case, subtract the hole contribution within R to the angularity
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['angularity_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['angularity_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['angularity_alice']['pt'][-1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['angularity_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['angularity_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:
                            for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                                kappa=1
                                if jet_collection_label in ['', '_shower_recoil', '_constituent_subtraction']:
                                    lambda_alpha = fjext.lambda_beta_kappa(jet, alpha, kappa, jetR)
                                elif jet_collection_label in ['_negative_recombiner']:
                                    lambda_alpha = 0
                                    for hadron in jet.constituents():
                                        if hadron.user_index() > 0:
                                            lambda_alpha += hadron.pt() / jet_pt * np.power(hadron.delta_R(jet)/jetR, alpha)
                                if jet_collection_label in ['_shower_recoil']:
                                    self.observable_dict_event[f'inclusive_chjet_angularity_alice_R{jetR}_alpha{alpha}{jet_collection_label}_unsubtracted'].append([jet_pt, lambda_alpha])
                                if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                    for hadron in holes_in_jet:
                                        if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                            continue
                                        lambda_alpha -= hadron.pt() / jet_pt  * np.power(hadron.delta_R(jet)/jetR, alpha)
                                self.observable_dict_event[f'inclusive_chjet_angularity_alice_R{jetR}_alpha{alpha}{jet_collection_label}'].append([jet_pt, lambda_alpha])

        if self.sqrts == 2760:

            # ALICE charged jet RAA
            if self.centrality_accepted(self.inclusive_chjet_observables['pt_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['pt_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['pt_alice']['pt'][1]
                if jetR in self.inclusive_chjet_observables['pt_alice']['jet_R']:
                    if abs(jet.eta()) < (self.inclusive_chjet_observables['pt_alice']['eta_cut']):
                        if pt_min < jet_pt < pt_max:

                            # Check leading track requirement
                            accept_jet = False
                            for constituent in jet.constituents():
                                if constituent.pt() > self.inclusive_chjet_observables['pt_alice']['leading_track_min_pt']:
                                    # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                    if abs(pid_hadrons_positive[np.abs(constituent.user_index())-1]) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                        accept_jet = True
                            if accept_jet:
                                self.observable_dict_event[f'inclusive_chjet_pt_alice_R{jetR}{jet_collection_label}'].append(jet_pt)
                                if jet_collection_label in ['_shower_recoil']:
                                    self.observable_dict_event[f'inclusive_chjet_pt_alice_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

            # g
            #   Hole treatment:
            #    - For show_recoil case, subtract the hole contribution within R to the angularity (also store unsubtracted case)
            #    - For negative_recombiner case, subtract the hole contribution within R to the angularity
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['g_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['g_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['g_alice']['pt'][1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['g_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['g_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:
                            g = 0
                            for constituent in jet.constituents():
                                if jet_collection_label in ['_negative_recombiner'] and constituent.user_index() < 0:
                                    continue
                                g += constituent.pt() / jet_pt * constituent.delta_R(jet)
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_chjet_g_alice_R{jetR}{jet_collection_label}_unsubtracted'].append(g)
                            if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                for hadron in holes_in_jet:
                                    if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                        continue
                                    g -= hadron.pt() / jet_pt * hadron.delta_R(jet)
                            self.observable_dict_event[f'inclusive_chjet_g_alice_R{jetR}{jet_collection_label}'].append(g)

            # pTD
            #   Hole treatment:
            #    - For show_recoil case, subtract the hole contribution within R to the angularity (also store unsubtracted case)
            #    - For negative_recombiner case, subtract the hole contribution within R to the angularity
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['ptd_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['ptd_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['ptd_alice']['pt'][1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['ptd_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['ptd_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:
                            sum = 0
                            for constituent in jet.constituents():
                                if jet_collection_label in ['_negative_recombiner'] and constituent.user_index() < 0:
                                    continue
                                sum += np.power(constituent.pt(), 2)
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_chjet_ptd_alice_R{jetR}{jet_collection_label}_unsubtracted'].append(np.sqrt(sum) / jet_pt)
                            if jet_collection_label in ['_shower_recoil', '_negative_recombiner']:
                                for hadron in holes_in_jet:
                                    if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() > 0 :
                                        continue
                                    sum -= np.power(hadron.pt(), 2)
                            self.observable_dict_event[f'inclusive_chjet_ptd_alice_R{jetR}{jet_collection_label}'].append(np.sqrt(sum) / jet_pt)

            # Jet mass
            #   Hole treatment:
            #    - For show_recoil case, subtract recoils within R from four-vector (also store unsubtracted case)
            #    - For negative_recombiner case, no subtraction is needed
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['mass_alice']['centrality']):
                pt_min = self.inclusive_chjet_observables['mass_alice']['pt'][0]
                pt_max = self.inclusive_chjet_observables['mass_alice']['pt'][-1]
                if abs(jet.eta()) < (self.inclusive_chjet_observables['mass_alice']['eta_cut_R'] - jetR):
                    if jetR in self.inclusive_chjet_observables['mass_alice']['jet_R']:
                        if pt_min < jet_pt < pt_max:
                            jet_mass = jet.m()
                            if jet_collection_label in ['_shower_recoil']:
                                self.observable_dict_event[f'inclusive_chjet_mass_alice_R{jetR}{jet_collection_label}_unsubtracted'].append([jet_pt, jet_mass])
                                # Subtract hole four vectors from the original jet, and then take the mass
                                jet_for_mass_calculation = fj.PseudoJet(jet)    # Avoid modifying the original jet.
                                for hadron in holes_in_jet:
                                    jet_for_mass_calculation -= hadron
                                # NOTE: Since we haven't assigned to `jet_mass` yet, it still contains the unsubtracted mass
                                jet_mass = jet_for_mass_calculation.m()
                            self.observable_dict_event[f'inclusive_chjet_mass_alice_R{jetR}{jet_collection_label}'].append([jet_pt, jet_mass])

        elif self.sqrts == 200:

            # STAR RAA
            if self.centrality_accepted(self.inclusive_chjet_observables['pt_star']['centrality']):
                pt_min = self.inclusive_chjet_observables['pt_star']['pt'][0]
                pt_max = 100. # Open upper bound
                if jetR in self.inclusive_chjet_observables['pt_star']['jet_R']:
                    if abs(jet.eta()) < (self.inclusive_chjet_observables['pt_star']['eta_cut_R'] - jetR):
                        if pt_min < jet_pt < pt_max:

                            # Check leading track requirement
                            min_leading_track_pt = 5.

                            accept_jet = False
                            for constituent in jet.constituents():
                                if constituent.pt() > min_leading_track_pt:
                                    # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                                    if abs(pid_hadrons_positive[np.abs(constituent.user_index())-1]) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                                        accept_jet = True
                            if accept_jet:
                                self.observable_dict_event[f'inclusive_chjet_pt_star_R{jetR}{jet_collection_label}'].append(jet_pt)
                                if jet_collection_label in ['_shower_recoil']:
                                    self.observable_dict_event[f'inclusive_chjet_pt_star_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

    # ---------------------------------------------------------------
    # Fill inclusive full jet observables
    # ---------------------------------------------------------------
    def fill_charged_jet_groomed_observables(self, grooming_setting, jet, jet_pt, jetR, jet_collection_label=''):

        # Construct groomed jet

        # For negative_recombiner case, we set the negative recombiner also for the C/A reclustering
        jet_def = fj.JetDefinition(fj.cambridge_algorithm, jetR)
        if jet_collection_label in ['_negative_recombiner']:
            recombiner = fjext.NegativeEnergyRecombiner()
            jet_def.set_recombiner(recombiner)
        gshop = fjcontrib.GroomerShop(jet, jet_def)

        zcut = grooming_setting['zcut']
        beta = grooming_setting['beta']
        jet_groomed_lund = gshop.soft_drop(beta, zcut, jetR)

        # ALICE hardest kt
        #   Hole treatment:
        #    - For show_recoil case, correct the pt only
        #    - For negative_recombiner case, no subtraction is needed
        #    - For constituent_subtraction, no subtraction is needed
        # For DyG, we need to record regardless of whether it passes SD, so we look at that observable first,
        # and then proceed with the rest afterwards.
        if self.sqrts == 5020:
            if self.centrality_accepted(self.inclusive_chjet_observables['kt_alice']['centrality']):
                # We put both DyG and SD after this setting to avoid filling DyG multiple times
                # (since DyG isn't included in the set of grooming settings)
                if grooming_setting in self.inclusive_chjet_observables['kt_alice']['SoftDrop']:
                    pt_min = self.inclusive_chjet_observables['kt_alice']['pt'][0]
                    pt_max = self.inclusive_chjet_observables['kt_alice']['pt'][-1]
                    if abs(jet.eta()) < (self.inclusive_chjet_observables['kt_alice']['eta_cut_R'] - jetR):
                        if jetR in self.inclusive_chjet_observables['kt_alice']['jet_R']:
                            if pt_min < jet_pt < pt_max:
                                for a in self.inclusive_chjet_observables['kt_alice']['dynamical_grooming_a']:
                                    jet_dyg_lund = gshop.dynamical(a)
                                    ktg = jet_dyg_lund.kt()
                                    self.observable_dict_event[f'inclusive_chjet_kt_alice_R{jetR}_a{a}{jet_collection_label}'].append([jet_pt, ktg])

                                # Only fill if SD identified a splitting
                                if jet_groomed_lund:
                                    ktg = jet_groomed_lund.kt()
                                    self.observable_dict_event[f'inclusive_chjet_kt_alice_R{jetR}_zcut{zcut}_beta{beta}{jet_collection_label}'].append([jet_pt, ktg])

        if not jet_groomed_lund:
            return

        if self.sqrts == 5020:

            # Soft Drop zg and theta_g
            #   Hole treatment:
            #    - For show_recoil case, correct the pt only
            #    - For negative_recombiner case, no subtraction is needed
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['zg_alice']['centrality']):
                if grooming_setting in self.inclusive_chjet_observables['zg_alice']['SoftDrop']:
                    pt_min = self.inclusive_chjet_observables['zg_alice']['pt'][0]
                    pt_max = self.inclusive_chjet_observables['zg_alice']['pt'][1]
                    if abs(jet.eta()) < (self.inclusive_chjet_observables['zg_alice']['eta_cut_R'] - jetR):
                        if jetR in self.inclusive_chjet_observables['zg_alice']['jet_R']:
                            if pt_min < jet_pt < pt_max:
                                theta_g = jet_groomed_lund.Delta() / jetR
                                zg = jet_groomed_lund.z()
                                # Note: untagged jets will return negative value
                                self.observable_dict_event[f'inclusive_chjet_zg_alice_R{jetR}_zcut{zcut}_beta{beta}{jet_collection_label}'].append(zg)
                                self.observable_dict_event[f'inclusive_chjet_tg_alice_R{jetR}_zcut{zcut}_beta{beta}{jet_collection_label}'].append(theta_g)

            # ALICE groomed angularity
            #   Hole treatment:
            #    - For show_recoil case, correct the pt only
            #    - For negative_recombiner case, correct the pt only
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.inclusive_chjet_observables['angularity_alice']['centrality']):
                if grooming_setting in self.inclusive_chjet_observables['angularity_alice']['SoftDrop']:
                    pt_min = self.inclusive_chjet_observables['angularity_alice']['pt'][0]
                    pt_max = self.inclusive_chjet_observables['angularity_alice']['pt'][-1]
                    if abs(jet.eta()) < (self.inclusive_chjet_observables['angularity_alice']['eta_cut_R'] - jetR):
                        if jetR in self.inclusive_chjet_observables['angularity_alice']['jet_R']:
                            if pt_min < jet_pt < pt_max:
                                for alpha in self.inclusive_chjet_observables['angularity_alice']['alpha']:
                                    kappa=1
                                    lambda_alpha = fjext.lambda_beta_kappa(jet_groomed_lund.pair(), alpha, kappa, jetR)
                                    self.observable_dict_event[f'inclusive_chjet_angularity_alice_R{jetR}_alpha{alpha}_zcut{zcut}_beta{beta}{jet_collection_label}'].append([jet_pt, lambda_alpha])

    # ---------------------------------------------------------------
    # Fill semi-inclusive charged jet observables
    # ---------------------------------------------------------------
    def fill_semi_inclusive_chjet_observables(self, jets_selected, hadrons_for_jet_finding, hadrons_negative, jetR, jet_collection_label=''):

        if self.sqrts == 2760:

            # Define trigger classes for both traditional h-jet analysis and Nsubjettiness analysis
            hjet_low_trigger_range = self.semi_inclusive_chjet_observables['IAA_alice']['low_trigger_range']
            hjet_high_trigger_range = self.semi_inclusive_chjet_observables['IAA_alice']['high_trigger_range']
            nsubjettiness_low_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['low_trigger_range']
            nsubjettiness_high_trigger_range = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['high_trigger_range']

            pt_IAA = self.semi_inclusive_chjet_observables['IAA_alice']['pt']
            pt_dphi = self.semi_inclusive_chjet_observables['dphi_alice']['pt']
            pt_nsubjettiness = self.semi_inclusive_chjet_observables['nsubjettiness_alice']['pt']

            # Define Nsubjettiness calculators
            axis_definition = fjcontrib.KT_Axes()
            measure_definition = fjcontrib.UnnormalizedMeasure(1)
            n_subjettiness_calculator1 = fjcontrib.Nsubjettiness(1, axis_definition, measure_definition)
            n_subjettiness_calculator2 = fjcontrib.Nsubjettiness(2, axis_definition, measure_definition)

            for hadron in hadrons_for_jet_finding:

                if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() < 0:
                    continue

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
                        # NOTE: This will record the hadron trigger even if it's not used in the IAA. However,
                        #       this is fine because we account for the difference in low and high trigger ranges
                        #       when we construct the histograms.
                        if jetR == min(self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']):
                            self.observable_dict_event[f'semi_inclusive_chjet_alice_trigger_pt{jet_collection_label}'].append(hadron.pt())

                        # Search for recoil jets
                        for jet in jets_selected:
                            if abs(jet.eta()) < (self.semi_inclusive_chjet_observables['IAA_alice']['eta_cut_R'] - jetR):

                                if jet_collection_label in ['_shower_recoil']:
                                    # Get the corrected jet pt: shower+recoil-holes
                                    jet_pt_unsubtracted = jet.pt()
                                    jet_pt_holes = 0
                                    for temp_hadron in hadrons_negative:
                                        if jet.delta_R(temp_hadron) < jetR:
                                            jet_pt_holes += temp_hadron.pt()
                                    jet_pt = jet_pt_unsubtracted - jet_pt_holes
                                else:
                                    jet_pt = jet_pt_unsubtracted = jet.pt()

                                # Jet yield and Delta phi
                                #   Hole treatment:
                                #    - For show_recoil case, correct the pt only (and also store unsubtracted pt)
                                #    - For negative_recombiner case, no subtraction is needed
                                #    - For constituent_subtraction, no subtraction is needed
                                if self.centrality_accepted(self.semi_inclusive_chjet_observables['IAA_alice']['centrality']):
                                    if hjet_found_low:
                                        if jetR in self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']:
                                            if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                                if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                    self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_lowTrigger{jet_collection_label}'].append(jet_pt)
                                                    if jet_collection_label in ['_shower_recoil']:
                                                        self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_lowTrigger{jet_collection_label}_unsubtracted'].append(jet_pt_unsubtracted)

                                        if jetR in self.semi_inclusive_chjet_observables['dphi_alice']['jet_R']:
                                            if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_dphi_alice_R{jetR}_lowTrigger{jet_collection_label}'].append(np.abs(hadron.delta_phi_to(jet)))

                                    if hjet_found_high:
                                        if jetR in self.semi_inclusive_chjet_observables['IAA_alice']['jet_R']:
                                            if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                                if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                    self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_highTrigger{jet_collection_label}'].append(jet_pt)
                                                    if jet_collection_label in ['_shower_recoil']:
                                                        self.observable_dict_event[f'semi_inclusive_chjet_IAA_alice_R{jetR}_highTrigger{jet_collection_label}_unsubtracted'].append(jet_pt_unsubtracted)

                                        if jetR in self.semi_inclusive_chjet_observables['dphi_alice']['jet_R']:
                                            if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_dphi_alice_R{jetR}_highTrigger{jet_collection_label}'].append(np.abs(hadron.delta_phi_to(jet)))

                                # Nsubjettiness
                                #   Hole treatment:
                                #    - For show_recoil case, correct the pt only
                                #    - For negative_recombiner case, no subtraction is needed
                                #    - For constituent_subtraction, no subtraction is needed
                                if self.centrality_accepted(self.semi_inclusive_chjet_observables['nsubjettiness_alice']['centrality']):
                                    if jetR in self.semi_inclusive_chjet_observables['nsubjettiness_alice']['jet_R']:
                                        if nsubjettiness_found_low:
                                            if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                                if pt_nsubjettiness[0] < jet_pt < pt_nsubjettiness[1]:
                                                    # We use the unsubtracted jet pt here, since the Nsubjettiness is calculated
                                                    # including recoils (but without hole subtraction) in the reclustering.
                                                    # Not ideal but not sure of an immediate better solution.
                                                    tau1 = n_subjettiness_calculator1.result(jet)/jet_pt_unsubtracted
                                                    tau2 = n_subjettiness_calculator2.result(jet)/jet_pt_unsubtracted
                                                    if tau1 > 1e-3:
                                                        self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_alice_R{jetR}_lowTrigger{jet_collection_label}'].append(tau2/tau1)

                                        if nsubjettiness_found_high:
                                            if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                                if pt_nsubjettiness[0] < jet_pt < pt_nsubjettiness[1]:
                                                    tau1 = n_subjettiness_calculator1.result(jet)/jet_pt_unsubtracted
                                                    tau2 = n_subjettiness_calculator2.result(jet)/jet_pt_unsubtracted
                                                    if tau1 > 1e-6:
                                                        self.observable_dict_event[f'semi_inclusive_chjet_nsubjettiness_alice_R{jetR}_highTrigger{jet_collection_label}'].append(tau2/tau1)

        if self.sqrts == 200:

            hjet_trigger_range = self.semi_inclusive_chjet_observables['IAA_star']['trigger_range']
            pt_IAA = self.semi_inclusive_chjet_observables['IAA_star']['pt']
            pt_dphi = self.semi_inclusive_chjet_observables['dphi_star']['pt']

            for hadron in hadrons_for_jet_finding:

                if jet_collection_label in ['_negative_recombiner'] and hadron.user_index() < 0:
                    continue

                if abs(hadron.eta()) < self.semi_inclusive_chjet_observables['IAA_star']['hadron_eta_cut']:

                    # Search for hadron trigger
                    found_trigger = False
                    if hjet_trigger_range[0] < hadron.pt() < hjet_trigger_range[1]:
                        found_trigger = True

                    if found_trigger:

                        # Record hadron pt for trigger normalization
                        if jetR == min(self.semi_inclusive_chjet_observables['IAA_star']['jet_R']):
                            self.observable_dict_event[f'semi_inclusive_chjet_star_trigger_pt{jet_collection_label}'].append(hadron.pt())

                        # Search for recoil jets
                        for jet in jets_selected:
                            if abs(jet.eta()) < (self.semi_inclusive_chjet_observables['IAA_star']['eta_cut_R'] - jetR):

                                if jet_collection_label in ['_shower_recoil']:
                                    # Get the corrected jet pt: shower+recoil-holes
                                    jet_pt_unsubtracted = jet.pt()
                                    jet_pt_holes = 0
                                    for temp_hadron in hadrons_negative:
                                        if jet.delta_R(temp_hadron) < jetR:
                                            jet_pt_holes += temp_hadron.pt()
                                    jet_pt = jet_pt_unsubtracted - jet_pt_holes
                                else:
                                    jet_pt = jet_pt_unsubtracted = jet.pt()

                                # Jet yield and Delta phi
                                if self.centrality_accepted(self.semi_inclusive_chjet_observables['IAA_star']['centrality']):

                                    if jetR in self.semi_inclusive_chjet_observables['IAA_star']['jet_R']:
                                            if np.abs(jet.delta_phi_to(hadron)) > (np.pi - 0.6):
                                                if pt_IAA[0] < jet_pt < pt_IAA[1]:
                                                    self.observable_dict_event[f'semi_inclusive_chjet_IAA_star_R{jetR}{jet_collection_label}'].append(jet_pt)
                                                    if jet_collection_label in ['_shower_recoil']:
                                                        self.observable_dict_event[f'semi_inclusive_chjet_IAA_star_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_unsubtracted)

                                    if jetR in self.semi_inclusive_chjet_observables['dphi_star']['jet_R']:
                                            if pt_dphi[0] < jet_pt < pt_dphi[1]:
                                                self.observable_dict_event[f'semi_inclusive_chjet_dphi_star_R{jetR}{jet_collection_label}'].append(np.abs(hadron.delta_phi_to(jet)))

    # ---------------------------------------------------------------
    # Fill dijet observables
    # ---------------------------------------------------------------
    def fill_dijet_observables(self, jets_selected, fj_hadrons_negative, jetR, jet_collection_label=''):

        if self.sqrts == 2760:

            # ATLAS xj
            #   Hole treatment:
            #    - For show_recoil case, correct jet pt by subtracting holes within R
            #    - For negative_recombiner case, no subtraction is needed
            #    - For constituent_subtraction, no subtraction is needed
            if self.centrality_accepted(self.dijet_observables['xj_atlas']['centrality']):
                if jetR in self.dijet_observables['xj_atlas']['jet_R']:

                    # First, find jets passing kinematic cuts
                    jet_candidates = []
                    for jet in jets_selected:

                        if jet_collection_label in ['_shower_recoil']:
                            # Get the corrected jet pt by subtracting the negative recoils within R
                            jet_pt = jet.pt()
                            if fj_hadrons_negative:
                                for temp_hadron in fj_hadrons_negative:
                                    if jet.delta_R(temp_hadron) < jetR:
                                        jet_pt -= temp_hadron.pt()
                        else:
                             jet_pt = jet.pt()

                        if jet_pt > self.dijet_observables['xj_atlas']['pt_subleading_min']:
                            if np.abs(jet.eta()) < self.dijet_observables['xj_atlas']['eta_cut']:
                                jet_candidates.append(jet)

                    # Find the leading two jets
                    leading_jet, leading_jet_pt, i_leading_jet = self.leading_jet(jet_candidates, fj_hadrons_negative, jetR)
                    if leading_jet:
                        del jet_candidates[i_leading_jet]
                        subleading_jet, subleading_jet_pt, _ = self.leading_jet(jet_candidates, fj_hadrons_negative, jetR)
                        if subleading_jet:
                            if np.abs(leading_jet.delta_phi_to(subleading_jet)) > 7*np.pi/8:
                                pt_min = self.dijet_observables['xj_atlas']['pt'][0]
                                if leading_jet_pt > pt_min:
                                    xj = subleading_jet_pt / leading_jet_pt
                                    self.observable_dict_event[f'dijet_xj_atlas_R{jetR}{jet_collection_label}'].append([leading_jet_pt, xj])

    #---------------------------------------------------------------
    # Return leading jet (or subjet)
    #---------------------------------------------------------------
    def leading_jet(self, jets, fj_hadrons_negative, jetR):

        leading_jet = None
        leading_jet_pt = 0.
        i_leading = 0
        for i,jet in enumerate(jets):

            # Get the corrected jet pt by subtracting the negative recoils within R
            jet_pt = jet.pt()
            for temp_hadron in fj_hadrons_negative:
                if jet.delta_R(temp_hadron) < jetR:
                    jet_pt -= temp_hadron.pt()

            if not leading_jet:
                leading_jet = jet
                leading_jet_pt = jet_pt
                i_leading = i

            if jet_pt > leading_jet_pt:
                leading_jet = jet
                leading_jet_pt = jet_pt
                i_leading = i

        return leading_jet, leading_jet_pt, i_leading

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
