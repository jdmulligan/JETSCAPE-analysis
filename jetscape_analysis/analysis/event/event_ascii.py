#!/usr/bin/env python3

"""
  Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
from event import event_base

################################################################
class event_ascii(event_base.event_base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, event="", **kwargs):
        super(event_ascii, self).__init__(**kwargs)

        self.event = event
