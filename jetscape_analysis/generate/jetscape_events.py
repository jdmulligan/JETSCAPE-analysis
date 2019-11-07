#!/usr/bin/env python3

"""
  Class to launch the generation of JETSCAPE events over a set of pt-hat bins.
  (You can easily adapt this to loop over any set of JETSCAPE parameter values).

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

# Base class
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

        self.pt_hat_bins = config["pt_hat_bins"]

        self.xml_user_file = config["xml_user_file"]
        self.xml_master_file = config["xml_master_file"]

    # ---------------------------------------------------------------
    # Main processing function
    # ---------------------------------------------------------------
    def generate_jetscape_events(self):

        # Loop through pT-hat bins and create directory structure with XML files for each bin
        for bin, pt_hat_min in enumerate(self.pt_hat_bins):

            # Set min,max of pT-hat bin
            if bin < (len(self.pt_hat_bins) - 1):
                pt_hat_max = self.pt_hat_bins[bin + 1]
            else:
                continue

            # Create outputDir for each bin
            output_dir_bin = "{}{}".format(self.output_dir, bin)
            if not output_dir_bin.endswith("/"):
                output_dir_bin = output_dir_bin + "/"
            if not os.path.exists(output_dir_bin):
                os.makedirs(output_dir_bin)

            # Set pT-hat values in Jetscape User XML configuration
            for line in fileinput.input(self.xml_user_file, inplace=True):
                if "pTHatMin" in line:
                    print("      <pTHatMin>{}</pTHatMin>".format(pt_hat_min))
                elif "pTHatMax" in line:
                    print("      <pTHatMax>{}</pTHatMax>".format(pt_hat_max))
                else:
                    print(line, end="")
            shutil.copyfile(self.xml_user_file, "{}{}".format(output_dir_bin, "jetscape_user.xml"))
            shutil.copyfile(
                self.xml_master_file, "{}{}".format(output_dir_bin, "jetscape_master.xml"),
            )

        # Loop through pt-hat bins and call Jetscape executable, and write output to pT-hat bin directory
        for bin, pt_hat_min in enumerate(self.pt_hat_bins):

            # Only iterate over lower bin edges
            if bin < (len(self.pt_hat_bins) - 1):
                pt_hat_max = self.pt_hat_bins[bin + 1]
                print("Generating pt-hat: {} - {} ...".format(pt_hat_min, pt_hat_max))
            else:
                continue

            # cd into bin directory in order to write Jetscape output there
            output_dir_bin = "{}{}".format(self.output_dir, bin)
            os.chdir(output_dir_bin)

            # Run Jetscape executable
            logfile_name = os.path.join(output_dir_bin, "log_{}.txt".format(bin))
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
