#!/usr/bin/env python3

""" Skim a large ascii file and split into smaller files
"""

from __future__ import print_function

# General
import os
import sys
import argparse
sys.path.append('.')

from jetscape_analysis.analysis.reader import parse_ascii
from jetscape_analysis.base import common_base

################################################################
class SkimAscii(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_file="", output_dir="", **kwargs):
        super(SkimAscii, self).__init__(**kwargs)
        self.input_file = input_file
        self.output_dir = output_dir
        
        self.event_id = 0
        self.events_per_chunk = 5000
    
    # ---------------------------------------------------------------
    # Main processing function for a single pt-hat bin
    # ---------------------------------------------------------------
    def skim(self):

        # Create reader class for each chunk of events, and iterate through each chunk
        # The parser returns an awkward array of events
        parse_ascii.parse_to_parquet(base_output_filename=self.output_dir,
                                     store_only_necessary_columns=True,
                                     input_filename=self.input_file,
                                     events_per_chunk=self.events_per_chunk)

##################################################################
if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser(description="Generate JETSCAPE events")
    parser.add_argument(
        "-i",
        "--inputFile",
        action="store",
        type=str,
        metavar="inputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/test.out",
        help="Input directory containing JETSCAPE output files",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        type=str,
        metavar="outputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Output directory and filename template for output to be written to",
    )

    # Parse the arguments
    args = parser.parse_args()

    # If invalid inputDir is given, exit
    if not os.path.exists(args.inputFile):
        print('File "{0}" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    analysis = SkimAscii(input_file=args.inputFile, output_dir=args.outputDir)
    analysis.skim()
