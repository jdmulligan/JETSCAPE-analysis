"""
Class to analyze a single JETSCAPE parquet output file,
and write out a new parquet file containing jets and jet constituents

Adapted from analyze_events_STAT.py

.. code-author:: Raymond Ehlers <raymond.ehlers@cern.ch>, LBL/UCB
"""

from __future__ import annotations

import argparse

# General
from collections import defaultdict
from pathlib import Path

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import numpy as np
import yaml
import awkward as ak

from jetscape_analysis.analysis import analyze_events_base_jet_constituents


class AnalyzeJetscapeEvents_Constituents(analyze_events_base_jet_constituents.AnalyzeJetscapeEvents_BaseJetConstituents):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super().__init__(config_file=config_file,
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
        with Path(self.config_file).open() as stream:
            config = yaml.safe_load(stream)

        self.sqrts = config['sqrt_s']
        self.output_file = 'observables'
        # Update the output_file to contain the labeling in the final_state_hadrons file.
        # We use this naming convention as the flag for whether we should attempt to rename it.
        if "final_state_hadrons" in self.input_file_hadrons:
            _input_filename = Path(self.input_file_hadrons).name
            # The filename will be something like "observables_0000_00.parquet", assuming
            # that the original name was "observables"
            self.output_file = _input_filename.replace("final_state_hadrons", self.output_file)
            #print(f'Updated output_file name to "{self.output_file}" in order to add identifying indices.')

        # Load observable blocks
        self.inclusive_jet_observables = config["inclusive_jet"]

        # General jet finding parameters
        self.jet_R = config['jet_R']
        self.min_jet_pt = config['min_jet_pt']
        self.max_jet_y = config['max_jet_y']

        # MG: Setup whatever general settings you might need here

        # If AA, set different options for hole subtraction treatment
        if self.is_AA:
            self.jet_collection_labels = config['jet_collection_labels']
        else:
            self.jet_collection_labels = ['']

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):
        # MG: Start here for analyzing an event. !

        # Initialize a dictionary that will store a list of calculated values for each output observable
        # MG: You can store the jet constituents here. Try something similar to do what you did the earlier pythia studies
        self.observable_dict_event = defaultdict(list)

        # MG: You don't need to worry too much about how this works to get started. Just jump down to `fill_jet_observables`.
        # Create list of fastjet::PseudoJets (separately for jet shower particles and holes)
        fj_hadrons_positive, pid_hadrons_positive = self.fill_fastjet_constituents(event, select_status='+')
        fj_hadrons_negative, pid_hadrons_negative = self.fill_fastjet_constituents(event, select_status='-')

        # Create list of charged particles
        fj_hadrons_positive_charged, pid_hadrons_positive_charged = self.fill_fastjet_constituents(event, select_status='+',
                                                                     select_charged=True)
        fj_hadrons_negative_charged, pid_hadrons_negative_charged = self.fill_fastjet_constituents(event, select_status='-',
                                                                     select_charged=True)

        # Fill jet observables
        # This loops over collections of jets, which is to say, different treatments of jets.
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

        # MG: We setup jet finding here, and then pass to `find_jets_and_fill`
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

            ## Charged jets
            #self.find_jets_and_fill(hadrons_for_jet_finding_charged, hadrons_negative_charged,
            #                        pid_hadrons_positive_charged, pid_hadrons_negative_charged,
            #                        jet_def, jet_selector, jetR, jet_collection_label, full_jet=False)

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


    # ---------------------------------------------------------------
    # Fill inclusive jet observables
    # ---------------------------------------------------------------
    def analyze_inclusive_jet(self, jet, hadrons_for_jet_finding, hadrons_negative,
                              pid_hadrons_positive, pid_hadrons_negative,
                              jetR, full_jet=True, jet_collection_label=''):

        # MG: This is where we actually start setting up to measure observables / extract jets.

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

    # ---------------------------------------------------------------
    # Fill inclusive full jet observables
    # ---------------------------------------------------------------
    def fill_full_jet_ungroomed_observables(self, jet, hadrons_for_jet_finding, holes_in_jet,
                                            pid_hadrons_positive, pid_hadrons_negative,
                                            jet_pt, jet_pt_uncorrected, jetR, jet_collection_label=''):
        # MG: This is just an example of how you could configure an analysis. You can do whatever you want!
        if self.centrality_accepted(self.inclusive_jet_observables['jet_constituents']['centrality']):
            # MG: Determine if the jet is accepted and store the constituents...
            if abs(jet.eta()) < (self.inclusive_jet_observables['jet_constituents']['eta_cut_R'] - jetR):

                # Store constituents
                constituents = {}
                constituents['px'] = []
                constituents['py'] = []
                constituents['pz'] = []
                constituents['E'] = []
                for constituent in jet.constituents():
                    constituents['px'].append(constituent.px())
                    constituents['py'].append(constituent.py())
                    constituents['pz'].append(constituent.pz())
                    constituents['E'].append(constituent.E())

                constituents = ak.zip(constituents)

                self.observable_dict_event[f'inclusive_jet_constituents_R{jetR}{jet_collection_label}'].append(constituents)



        ## ALICE RAA
        ##   Hole treatment:
        ##    - For RAA, all jet collections can be filled from the corrected jet pt
        ##    - In the shower_recoil case, we also fill the unsubtracted jet pt
        #if self.centrality_accepted(self.inclusive_jet_observables['pt_alice']['centrality']):
        #    pt_min, pt_max = self.inclusive_jet_observables['pt_alice']['pt']
        #    if jetR in self.inclusive_jet_observables['pt_alice']['jet_R']:
        #        if abs(jet.eta()) < (self.inclusive_jet_observables['pt_alice']['eta_cut_R'] - jetR):
        #            if pt_min < jet_pt < pt_max:

        #                # Check leading track requirement
        #                if jetR == 0.2:
        #                    min_leading_track_pt = 5.
        #                else:
        #                    min_leading_track_pt = 7.

        #                accept_jet = False
        #                acceptable_hadrons = [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]
        #                for constituent in jet.constituents():
        #                    if constituent.pt() > min_leading_track_pt:
        #                        # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
        #                        if abs(pid_hadrons_positive[np.abs(constituent.user_index())-1]) in acceptable_hadrons:
        #                            accept_jet = True
        #                if accept_jet:
        #                    self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}{jet_collection_label}'].append(jet_pt)
        #                    if jet_collection_label in ['_shower_recoil']:
        #                        self.observable_dict_event[f'inclusive_jet_pt_alice_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

        ## ATLAS RAA
        #if self.centrality_accepted(self.inclusive_jet_observables['pt_atlas']['centrality']):
        #    pt_min = self.inclusive_jet_observables['pt_atlas']['pt'][0]
        #    pt_max = self.inclusive_jet_observables['pt_atlas']['pt'][1]
        #    if jetR in self.inclusive_jet_observables['pt_atlas']['jet_R']:
        #        if abs(jet.rap()) < self.inclusive_jet_observables['pt_atlas']['y_cut']:
        #            if pt_min < jet_pt < pt_max:
        #                self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}{jet_collection_label}'].append(jet_pt)
        #                if jet_collection_label in ['_shower_recoil']:
        #                    self.observable_dict_event[f'inclusive_jet_pt_atlas_R{jetR}{jet_collection_label}_unsubtracted'].append(jet_pt_uncorrected)

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

            if fj_hadrons_negative:
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
            msg = f'failed to compute charge of pid {pid}'
            raise ValueError(msg)


def main() -> None:
    """Main entry point
    """
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
    if not Path(args.configFile).exists():
        msg = f'File "{args.configFile}" does not exist! Exiting!'
        raise ValueError(msg)

    # If invalid inputDir is given, exit
    if not Path(args.inputFile).exists():
        msg = f'File "{args.inputFile}" does not exist! Exiting!'
        raise ValueError(msg)

    analysis = AnalyzeJetscapeEvents_Constituents(
        config_file=args.configFile,
        input_file=args.inputFile,
        output_dir=args.outputDir
    )
    analysis.analyze_jetscape_events()


if __name__ == "__main__":
    main()