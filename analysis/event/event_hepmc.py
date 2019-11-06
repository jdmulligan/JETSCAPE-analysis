#!/usr/bin/env python3

"""
  HepMC Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
from event import event_base

################################################################
class event_hepmc(event_base.event_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, event='', **kwargs):
    super(event_hepmc, self).__init__(**kwargs)
  
    self.event = event

  #---------------------------------------------------------------
  # Get list of hadrons.
  # Final state hadrons (from jet + bulk) are stored as outgoing particles in a disjoint vertex with t = 100
  #---------------------------------------------------------------
  def get_hadrons(self):
    
    for vertex in self.event.vertices:
      
      vertex_time = vertex.position.t
      if abs(vertex_time - 100) < 1e-3:
        final_state_particles = vertex.particles_out
    
    # Remove neutrinos
    hadrons = []
    for particle in final_state_particles:
      
      pid = particle.pid
      if pid!=12 and pid!=14 and pid!=16:
        hadrons.append(particle)

    return hadrons
