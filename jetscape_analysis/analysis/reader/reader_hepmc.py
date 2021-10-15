#!/usr/bin/env python3

"""
  HepMC3 reader class, based on pyhepmc_ng
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# General
import pyhepmc_ng

# Event class
from event import event_hepmc

# Base class
from reader import reader_base

################################################################
class ReaderHepMC(reader_base.ReaderBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file="", **kwargs):
        super(ReaderHepMC, self).__init__(**kwargs)

        # Create reader
        self.reader = pyhepmc_ng.ReaderAscii(input_file)
        if self.reader.failed():
            raise ValueError("[error] unable to read from {}".format(input_file))

    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def next_event(self):

        event = pyhepmc_ng.GenEvent()
        self.reader.read_event(event)

        if self.reader.failed():
            return False
        else:
            ev = event_hepmc.EventHepMC(event)
            return ev
