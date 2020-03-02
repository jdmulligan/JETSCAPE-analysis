#!/usr/bin/env python3

"""
  Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
from event import event_base

################################################################
class EventAscii(event_base.EventBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, event="", **kwargs):
        super(EventAscii, self).__init__(**kwargs)

        self.event = event
