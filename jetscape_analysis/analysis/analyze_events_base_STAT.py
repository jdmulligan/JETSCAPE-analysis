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
import pickle

# Analysis
import itertools
import pandas as pd
import numpy as np
import ROOT

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
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        self.n_event_max = 100
            
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

    # ---------------------------------------------------------------
    # Analyze event chunk
    # ---------------------------------------------------------------
    def analyze_event_chunk(self, df_event_chunk):

        # Loop through events
        for i,event in df_event_chunk.iterrows():

            if i % 1000 == 0:
                print(f'event: {i}')
                
            if i > self.n_event_max:
                return
                
            # Store dictionary of all observables for the event
            self.observable_dict_event = {}
        
            # Fill event cross-section
            xsec = 0.1
            self.observable_dict_event['xsec'] = xsec

            # Call user-defined function to analyze event
            self.analyze_event(event)

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_output_objects(self):
    
        # Initialize list to store observables
        # Each entry in the list stores a dict for a given event
        self.output_event_list = []

    # ---------------------------------------------------------------
    # Save output event list into a dataframe
    # ---------------------------------------------------------------
    def write_output_objects(self):
        
        self.output_dataframe = pd.DataFrame(self.output_event_list)
        with open(self.output_file, 'wb') as f:
            pickle.dump(self.output_dataframe, f)
        print(self.output_dataframe.keys())
        print(self.output_dataframe)

    # ---------------------------------------------------------------
    # Fill hadrons into vector of fastjet pseudojets
    #
    # By default, select all particles
    # If select_status='+', select only positive status particles
    # If select_status='-', select only positive status particles
    # ---------------------------------------------------------------
    def fill_fastjet_constituents(self, event, select_status=None, select_charged=False):
    
        if select_status == '-':
            status_mask = (event['status'] < 0)
        elif select_status == '+':
            status_mask = (event['status'] > -1)
        else:
            # Picked a value to make an all true mask. We don't select anything
            status_mask = event["status"] > -1e6

        # Default to an all true mask
        charged_mask = np.ones(len(status_mask)) > 0
        if select_charged:
            # NOTE: This is super inefficient - we can seriously accelerate this with numba.
            #       But this apparently works for now, so we leave it as is for now.
            # Create an all false mask. We'll fill it in with the charged constituents
            charged_mask = np.ones(len(event['particle_ID'])) < 0
            for i, pid_value in enumerate(event['particle_ID']):
                # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
                if np.abs(pid_value) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
                    charged_mask[i] = True

        full_mask = status_mask & charged_mask
        px = event['px'][full_mask]
        py = event['py'][full_mask]
        pz = event['pz'][full_mask]
        e = event['E'][full_mask]
        pid = event['particle_ID'][full_mask]

        # Create a vector of fastjet::PseudoJets from arrays of px,py,pz,e
        fj_particles = fjext.vectorize_px_py_pz_e(px, py, pz, e)

        # Set pid as user_index
        for i,p in enumerate(fj_particles):
            fj_particles[i].set_user_index(int(pid[i]))

        return fj_particles

    # ---------------------------------------------------------------
    # This function is called once per event
    # You must implement this
    # ---------------------------------------------------------------
    def analyze_event(self, event):
        raise NotImplementedError('You must implement analyze_event()!')
