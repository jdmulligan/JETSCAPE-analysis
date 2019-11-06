#!/usr/bin/env python3

"""
  Reader base class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
import common_base

################################################################
class reader_base(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', **kwargs):
    super(reader_base, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Generator (in pythonic sense) to loop over all events
  #---------------------------------------------------------------
  def __call__(self, n_events):
    
    for _ in range(0, n_events):
      yield self.next_event()

