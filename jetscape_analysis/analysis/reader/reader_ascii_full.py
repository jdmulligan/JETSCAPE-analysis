#!/usr/bin/env python3

"""
  Ascii reader class -- not currently working
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

import numpy as np

# Event class
from event import event_ascii
from reader import string_tokenizer
import pyhepmc_ng

# Base class
from reader import reader_base

################################################################
class ReaderAscii(reader_base.ReaderBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file="", **kwargs):
        super(ReaderAscii, self).__init__(**kwargs)
        
        self.input_file = input_file
        with open(input_file, 'r') as f:
            self.file = f.readlines()

        self.string_tokenizer = string_tokenizer.StringTokenizer()
        
        self.current_event = 0
        self.line_index = 0
        self.current_shower = 1
        self.p_shower = None
        self.p_showers = []
        self.node_vec = []
        self.edge_vec = []
        self.hadrons = []
        
    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def next_event(self):

        # Construct an empty ascii event
        event = event_ascii.EventAscii()
        
        # Fill the event
        return self.read_event(event)
    
    # ---------------------------------------------------------------
    # Get next event
    # Return event if successful, False if unsuccessful
    # ---------------------------------------------------------------
    def read_event(self, event):
        
        token = ''
        node_zero_counter = 0
    
        if self.current_event > 0:
            self.node_vec = []
            self.edge_vec = []
            self.p_showers = []
            self.hadrons = []

        print('event: {}'.format(self.current_event))
        for i, line in enumerate(self.file):
            
            print('index: {}'.format(self.line_index))

            self.line_index += 1
            if i < self.line_index + 1:
                continue
                
            self.string_tokenizer.set(line)
    
            if self.string_tokenizer.is_comment_entry():
                continue
    
            if self.string_tokenizer.is_event_entry():
            
                new_event = int(self.string_tokenizer.next())
                
                if new_event > 0:
                    print('new event: {}'.format(new_event))
                    print('current: {}'.format(self.current_event))

                # If it is the event header for the current event, continue
                # Otherwise, return the event
                if self.current_event != new_event:
                    self.current_event += 1
                    return event
                
            # not an event header -- done?
            if not self.string_tokenizer.is_graph_entry():
                continue
            
            # node?
            if self.string_tokenizer.is_node_entry():
                '''
                # catch starting node
                if self.string_tokenizer.is_node_zero():
                    node_zero_counter += 1

                    if node_zero_counter > self.current_shower:
                        self.node_vec = []
                        self.edge_vec = []
                        #pShowers.push_back(make_shared<PartonShower>());
                        #pShower=pShowers.back();
                        self.current_shower += 1
                        
                    self.add_node(line)
                    continue
                '''
                continue
                        
            # edge?
            if self.string_tokenizer.is_edge_entry():
                #self.add_edge(line)
                continue

      
            # else: hadron entry
            print('hadron entry')
            self.add_hadron(event, line)

    # ---------------------------------------------------------------
    # Add a hadron to the hadron list, given its ascii line
    # Construct it as a HepMC::GenParticle, so that it can be treated the same as in the HepMC case
    # ---------------------------------------------------------------
    def add_hadron(self, event, line):
    
        self.string_tokenizer.set(line)
        vars = []
        while self.string_tokenizer.done():
            token = self.string_tokenizer.next()
            if token != 'H':
                vars.append(token)
                
        label = int(vars[1])
        id = int(vars[2])
        stat = int(vars[3])
        pt = float(vars[4])
        eta = float(vars[5])
        phi = float(vars[6])
        e = float(vars[7])
        
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        
        four_vector = pyhepmc_ng.FourVector(px, py, pz, e)
        particle = pyhepmc_ng.GenParticle(four_vector, id, status)
        event.hadrons.append(particle)

    #---------------------------------------------------------------
    # Add node
    #---------------------------------------------------------------
    def add_node(self, s):
  
        token = ''

        self.string_tokenizer.set(s)
    
        vS = []
        while not self.string_tokenizer.done():
    
            token = self.string_tokenizer.next()
            if token.compare('V') != 0:
                vS.append(token)
            
              
        '''
        vector<string> vS;
        
        while (!strT.done())
          {
            token = strT.next();
            if(token.compare("V") != 0)
              vS.push_back(token);
        }

        nodeVec.push_back(pShower->new_vertex(make_shared<Vertex>(stod(vS[1]),stod(vS[2]),stod(vS[3]),stod(vS[4]))));
        '''
