#!/usr/bin/env python3

"""
  Ascii reader class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Event class
from event import event_ascii

# Base class
from reader import reader_base

################################################################
class ReaderAscii(reader_base.ReaderBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file="", **kwargs):
        super(ReaderAscii, self).__init__(**kwargs)

        # Create reader

    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def next_event(self):

        return False
