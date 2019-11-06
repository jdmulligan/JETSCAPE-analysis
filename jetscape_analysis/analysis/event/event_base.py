#!/usr/bin/env python3

"""
  Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
import common_base

################################################################
class event_base(common_base.common_base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(event_base, self).__init__(**kwargs)
