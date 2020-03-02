#!/usr/bin/env python3

"""
  HepMC Event class
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# Base class
from jetscape_analysis.analysis.event import event_base

################################################################
class EventHepMC(event_base.EventBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, event="", **kwargs):
        super(EventHepMC, self).__init__(**kwargs)

        self.event = event

    # ---------------------------------------------------------------
    # Get list of hadrons.
    # Final state hadrons (from jet + bulk) are stored as outgoing particles in a disjoint vertex with t = 100
    # ---------------------------------------------------------------
    def hadrons(self, min_track_pt=0.):

        for vertex in self.event.vertices:

            vertex_time = vertex.position.t
            if abs(vertex_time - 100) < 1e-3:
                final_state_particles = vertex.particles_out

        # Remove neutrinos
        hadrons = []
        for particle in final_state_particles:

            pid = particle.pid
            pt = particle.momentum.pt()
            if pid != 12 and pid != 14 and pid != 16:
                if pt > min_track_pt:
                    hadrons.append(particle)

        return hadrons

    # ---------------------------------------------------------------
    # Get list of final-state partons.
    # ---------------------------------------------------------------
    def final_partons(self):

        partons = []
        n_vertices = len(self.event.vertices)
        n_particles = len(self.event.particles)
        for particle in self.event.particles:

            status = particle.status
            n_children = len(particle.children)

            parent_vertex = particle.production_vertex
            end_vertex = particle.end_vertex
            parent_vertex_time = parent_vertex.position.t

            is_parton = False
            if abs(parent_vertex_time - 100) > 1e-3:
                is_parton = True

            is_final_parton = False
            if is_parton and not end_vertex:
                is_final_parton = True

            # Alternately: It seems that status=0 means the particle is a parton
            # Not sure what the nonzero status codes mean (62, 83, 84)

            if is_final_parton:
                partons.append(particle)

        return partons
