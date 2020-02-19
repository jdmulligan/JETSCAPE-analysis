#!/usr/bin/env python3

"""
  Class to launch the generation of JETSCAPE events over a set of pt-hat bins,
  and optionally over any set of additional JETSCAPE parameter values.
  If additional parameter values are specified, all combinations of the
  specified parameter values will be generated.

  Run from inside the JETSCAPE docker container with:
    python generate_jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml -o /my/outputdir

.. codeauthor:: James Mulligan <james.mulligan@berkeley.edu>, UC Berkeley
"""

from __future__ import print_function

import argparse
import fileinput
import os
import shutil
import subprocess
import sys
import yaml
import itertools

# Base class
sys.path.append('../..')
from jetscape_analysis.base import common_base

################################################################
class generate_jetscape_events(common_base.common_base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file="", output_dir="", **kwargs):
        super(generate_jetscape_events, self).__init__(**kwargs)
        self.config_file = config_file
        self.output_dir = output_dir

        # Create output dir
        if not self.output_dir.endswith("/"):
            self.output_dir = self.output_dir + "/"
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.initialize_config()

        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, "r") as stream:
            config = yaml.safe_load(stream)

        self.debug_level = config["debug_level"]

        self.xml_user_file = config["xml_user_file"]
        self.xml_master_file = config["xml_master_file"]
                
        self.parameter_scan_dict = config['parameter_scan']
        self.pt_hat_bins = self.parameter_scan_dict['pt_hat_bins']['values']

    # ---------------------------------------------------------------
    # Main processing function
    # ---------------------------------------------------------------
    def generate_jetscape_events(self):

        # Store list of parameter labels and formatting
        parameter_labels = [self.parameter_scan_dict[key]['label'] for key in self.parameter_scan_dict]
        parameter_spacings = [self.parameter_scan_dict[key]['spacing'] for key in self.parameter_scan_dict]

        # Create list of all combinations of parameters
        parameter_values = [self.parameter_scan_dict[key]['values'] for key in self.parameter_scan_dict]
        parameter_combinations = list(itertools.product(*parameter_values))
        
        # Remove that last pt-hat bin edge
        n_combinations_per_pthat = int(len(parameter_combinations)/len(self.pt_hat_bins))
        parameter_combinations = parameter_combinations[:-n_combinations_per_pthat]

        # Loop through all parameter combinations
        for index, parameter_combination in enumerate(parameter_combinations):
        
            pt_hat_bin = int(index / n_combinations_per_pthat)
            if pt_hat_bin < len(self.pt_hat_bins) - 1:
                pt_hat_min = self.pt_hat_bins[pt_hat_bin]
                pt_hat_max = self.pt_hat_bins[pt_hat_bin + 1]
            else:
                continue
            if index % n_combinations_per_pthat == 0:
                print('Generating pt-hat: {} - {} ...'.format(pt_hat_min, pt_hat_max))

            # Create label for output directory
            dir_label = ''
            for index, value in enumerate(parameter_combination):
                if index == 0:
                    continue
                if index != 1:
                    dir_label += '_'
                dir_label += parameter_labels[index]
                dir_label += str(value)
            if len(parameter_combination) > 1:
                print('    Generating {}'.format(dir_label))
                
            # Create outputDir for each bin
            output_dir_bin = '{}{}/{}'.format(self.output_dir, pt_hat_bin, dir_label)
            if not output_dir_bin.endswith("/"):
                output_dir_bin = output_dir_bin + "/"
            if not os.path.exists(output_dir_bin):
                os.makedirs(output_dir_bin)

            # Copy XML files to pt-hat bin directory
            xml_master_file_copy = "{}{}".format(output_dir_bin, "jetscape_master.xml")
            cmd = "rsync {} {}".format(self.xml_master_file, xml_master_file_copy)
            os.system(cmd)

            xml_user_file_copy = "{}{}".format(output_dir_bin, "jetscape_user.xml")
            cmd = "rsync {} {}".format(self.xml_user_file, xml_user_file_copy)
            os.system(cmd)
                    
            # Set parameters in the Jetscape User XML configuration
            for index, value in enumerate(parameter_combination):
            
                parameter_label = parameter_labels[index]
                parameter_spacing = parameter_spacings[index]
            
                # Set pt-hat
                if parameter_label == 'pt_hat_bins':
                
                    for line in fileinput.input(xml_user_file_copy, inplace=True):
                
                        if 'pTHatMin' in line:
                            print('{}<pTHatMin>{}</pTHatMin>'.format(parameter_spacing, pt_hat_min))
                        elif 'pTHatMax' in line:
                            print('{}<pTHatMax>{}</pTHatMax>'.format(parameter_spacing, pt_hat_max))
                        else:
                            print(line, end='')
                      
                # Set other parameters
                else:
                
                    for line in fileinput.input(xml_user_file_copy, inplace=True):
                    
                        if parameter_label in line:
                            print('{}<{}>{}</{}>'.format(parameter_spacing, parameter_label, value, parameter_label))
                        else:
                            print(line, end='')

            # cd into bin directory in order to write Jetscape output there
            os.chdir(output_dir_bin)

            # Run Jetscape executable
            logfile_name = os.path.join(output_dir_bin, "log_{}_{}.txt".format(pt_hat_bin, dir_label))
            with open(logfile_name, "w") as logfile:
                cmd = "/home/jetscape-user/JETSCAPE/build/runJetscape jetscape_user.xml jetscape_master.xml"
                subprocess.run(cmd, check=True, shell=True, stdout=logfile)
            os.chdir(self.output_dir)

##################################################################
if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser(description="Generate JETSCAPE events")
    parser.add_argument(
        "-c",
        "--configFile",
        action="store",
        type=str,
        metavar="configFile",
        default="/home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml",
        help="Path of config file for analysis",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        type=str,
        metavar="outputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Output directory for output to be written to",
    )

    # Parse the arguments
    args = parser.parse_args()

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File "{0}" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = generate_jetscape_events(config_file=args.configFile, output_dir=args.outputDir)
    analysis.generate_jetscape_events()
