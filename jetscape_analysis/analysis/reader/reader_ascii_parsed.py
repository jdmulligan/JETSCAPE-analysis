#!/usr/bin/env python3

"""
  Ascii reader class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

import os
import sys
import numpy as np
import pandas as pd

# Event class
from event import event_ascii

# Base class
from reader import reader_base

################################################################
class ReaderAsciiParsed(reader_base.ReaderBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, event_chunk_hadrons=None, event_chunk_partons=None, **kwargs):
        super(ReaderAsciiParsed, self).__init__(**kwargs)
        
        self.event_chunk_hadrons = event_chunk_hadrons
        self.event_chunk_partons = event_chunk_partons
        
        self.current_event = 0
        self.n_events = len(self.event_chunk_hadrons)

        if self.event_chunk_partons:
            if len(self.event_chunk_hadrons) != len(self.event_chunk_partons):
                sys.exit('Final state partons has {} events, but partons has {}.'.format(len(self.event_list_hadrons), len(self.event_list_partons)))
        
    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def next_event(self):

        if self.current_event < self.n_events:
            self.current_event += 1
            event_hadrons = self.event_chunk_hadrons[self.current_event-1]
            if self.event_chunk_partons:
                event_partons = self.event_chunk_partons[self.current_event-1]
            else:
                event_partons = None
            return event_ascii.EventAscii(event_hadrons, event_partons)
        else:
            sys.exit('Current event {} greater than total n_events {}'.format(self.current_event,
                                                                              self.n_events))
