#!/usr/bin/env python3

""" Base class to analyze a JETSCAPE output file

You should create a user class that inherits from this one. See analyze_events_STAT.py for an example.

The outputdir should contain a JETSCAPE output file in parquet format

See README for pre-requisites.

.. codeauthor:: James Mulligan <james.mulligan@berkeley.edu>, UC Berkeley
"""

from __future__ import print_function

# General
import os
import yaml
import time
from pathlib import Path
from numba import jit

# Analysis
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
from pathlib import Path

# Fastjet via python (from external library heppy)
import fjext

from jetscape_analysis.base import common_base

################################################################
class AnalyzeJetscapeEvents_BaseSTAT(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file="", input_file="", output_dir="", **kwargs):
        super(AnalyzeJetscapeEvents_BaseSTAT, self).__init__(**kwargs)

        self.config_file = config_file
        self.input_file_hadrons = input_file
        self.output_dir = Path(output_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)

            if 'n_event_max' in config:
                self.n_event_max = config['n_event_max']
            else:
                self.n_event_max = -1

        # Check whether pp or AA
        if 'PbPb' in self.input_file_hadrons or 'AuAu' in self.input_file_hadrons:
            self.is_AA = True
        else:
            self.is_AA = False

        # If AA, get centrality bin
        if self.is_AA:
            _final_state_hadrons_path = Path(self.input_file_hadrons)
            # For an example filename of "jetscape_PbPb_Run0005_5020_0001_final_state_hadrons_00.parquet",
            # - the run number is index 2
            _run_number = _final_state_hadrons_path.stem.split("_")[2]
            # - the file index is at index 4 (in the example, it extracts `1` as an int)
            _file_index = int(_final_state_hadrons_path.name.split('_')[4])
            run_info_path = _final_state_hadrons_path.parent / f"{_run_number}_info.yaml"
            with open(run_info_path, 'r') as f:
                _run_info = yaml.safe_load(f)
                centrality_string = _run_info["index_to_hydro_event"][_file_index].split('/')[0].split('_')
                # index of 1 and 2 based on an example entry of "cent_00_01"
                self.centrality = [int(centrality_string[1]), int(centrality_string[2])]

    # ---------------------------------------------------------------
    # Main processing function
    # ---------------------------------------------------------------
    def analyze_jetscape_events(self):

        print('Analyzing events ...')

        # Initialize output objects
        self.initialize_output_objects()

        # Read chunk of events into a dataframe
        # Fields: particle_ID, status, E, px, py, pz
        df_event_chunk = pd.read_parquet(self.input_file_hadrons)
        if self.n_event_max < 0:
            self.n_event_max = df_event_chunk.shape[0]

        # Iterate through events
        self.analyze_event_chunk(df_event_chunk)

        # Write analysis task output to ROOT file
        self.write_output_objects()

        print('Done!')

    # ---------------------------------------------------------------
    # Analyze event chunk
    # ---------------------------------------------------------------
    def analyze_event_chunk(self, df_event_chunk):

        # Loop through events
        start = time.time()
        weight_sum = 0.
        for i,event in df_event_chunk.iterrows():

            if i % 1000 == 0:
                print(f'event: {i}    (time elapsed: {time.time() - start} s)')

            if i > self.n_event_max:
                break

            # Store dictionary of all observables for the event
            self.observable_dict_event = {}

            # Call user-defined function to analyze event
            self.analyze_event(event)

            # Fill the observables dict to a new entry in the event list
            event_weight = event['event_weight']
            weight_sum += event_weight
            if self.event_has_entries(self.observable_dict_event):

                # Fill event cross-section weight
                self.observable_dict_event['event_weight'] = event_weight
                self.observable_dict_event['pt_hat'] = event['pt_hat']

                self.output_event_list.append(self.observable_dict_event)

        # Get total cross-section (same for all events at this point), weight sum, and centrality
        self.cross_section_dict['cross_section'] = event['cross_section']
        self.cross_section_dict['cross_section_error'] = event['cross_section_error']
        self.cross_section_dict['n_events'] = self.n_event_max
        self.cross_section_dict['weight_sum'] = weight_sum
        if self.is_AA:
            self.cross_section_dict['centrality_min'] = self.centrality[0]
            self.cross_section_dict['centrality_max'] = self.centrality[1]

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_output_objects(self):

        # Initialize list to store observables
        # Each entry in the list stores a dict for a given event
        self.output_event_list = []

        # Store also the total cross-section (one number per file)
        self.cross_section_dict = {}

    # ---------------------------------------------------------------
    # Save output event list into a dataframe
    # ---------------------------------------------------------------
    def event_has_entries(self, event_dict):

        return bool([obs for obs in event_dict.values() if obs != []])

    # ---------------------------------------------------------------
    # Check if event centrality is within observable's centrality
    # ---------------------------------------------------------------
    def centrality_accepted(self, observable_centrality_list):

        # AA
        if self.is_AA:

            for observable_centrality in observable_centrality_list:
                if self.centrality[0] >= observable_centrality[0]:
                    if self.centrality[1] <= observable_centrality[1]:
                        return True
            return False

        # pp
        else:
            return True

    # ---------------------------------------------------------------
    # Save output event list into a dataframe
    # ---------------------------------------------------------------
    def write_output_objects(self):

        # Convert to pandas, and then arrow.
        self.output_dataframe = pd.DataFrame(self.output_event_list)
        #self.output_dataframe = ak.Array(self.output_event_list)
        table = pa.Table.from_pandas(self.output_dataframe)

        # Write to parquet
        # Determine the types for improved compression when writing
        # See writing to parquet in the final state hadrons parser for more info.
        float_types = [np.float32, np.float64]
        float_columns = list(self.output_dataframe.select_dtypes(include=float_types).keys())
        other_columns = list(self.output_dataframe.select_dtypes(exclude=float_types).keys())
        # NOTE: As of 27 April 2021, this doesn't really work right because too many columns
        #       are of the "object" type. We may need to revise the output format to optimize
        #       the output size.
        print(f"float_columns: {float_columns}")
        print(f"other_columns: {other_columns}")
        pq.write_table(
            table, self.output_dir / self.output_file, compression="zstd",
            use_dictionary=other_columns,
            use_byte_stream_split=float_columns,
        )

        # Write cross-section to separate file
        cross_section_dataframe = pd.DataFrame(self.cross_section_dict, index=[0])
        cross_section_table = pa.Table.from_pandas(cross_section_dataframe)
        filename = self.output_file.replace('observables', 'cross_section')
        pq.write_table(cross_section_table, self.output_dir / filename, compression="zstd")

    # ---------------------------------------------------------------
    # Fill hadrons into vector of fastjet pseudojets
    #
    # By default, select all particles
    # If select_status='+', select only positive status particles
    # If select_status='-', select only positive status particles
    # ---------------------------------------------------------------
    def fill_fastjet_constituents(self, event, select_status=None, select_charged=False):

        # Construct indices according to particle status
        if select_status == '-':
            status_mask = (event['status'] < 0)
        elif select_status == '+':
            status_mask = (event['status'] > -1)
        else:
            # Picked a value to make an all true mask. We don't select anything
            status_mask = event["status"] > -1e6

        # Construct indices according to charge
        charged_mask = get_charged_mask(event['particle_ID'], select_charged)

        full_mask = status_mask & charged_mask
        px = event['px'][full_mask]
        py = event['py'][full_mask]
        pz = event['pz'][full_mask]
        e = event['E'][full_mask]
        pid = event['particle_ID'][full_mask]

        # Create a vector of fastjet::PseudoJets from arrays of px,py,pz,e
        fj_particles = fjext.vectorize_px_py_pz_e(px, py, pz, e)

        # Set pid as user_index
        [fj_particles[i].set_user_index(int(pid[i])) for i,_ in enumerate(fj_particles)]

        return fj_particles

    # ---------------------------------------------------------------
    # This function is called once per event
    # You must implement this
    # ---------------------------------------------------------------
    def analyze_event(self, event):
        raise NotImplementedError('You must implement analyze_event()!')

# ---------------------------------------------------------------
# Construct charged particle mask
# ---------------------------------------------------------------
@jit(nopython=True)
def get_charged_mask(pid, select_charged: bool):
    """ Create mask for selected a set of charged particles based on PID.

    Note:
        This function assumes that the same set of charged particles are selected
        for all charged-particle jets (ie. ALICE and STAR). Although the charged
        particle selections for some of the hadron observables vary between experiments,
        this seems like

    Args:
        pid: PID values associated with the charged particles in an event.
        select_charged: If True, actually select charged particles. If False,
            just return an all True mask (if for full jets).

    Returns:
        Mask selecting the particles.
    """
    # Default to an all true mask
    charged_mask = np.ones(len(pid)) > 0
    if select_charged:
        # Create an all false mask. We'll fill it in with the charged constituents
        charged_mask = np.ones(len(pid)) < 0
        for i, pid_value in enumerate(pid):
            # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if np.abs(pid_value) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                charged_mask[i] = True

    return charged_mask
