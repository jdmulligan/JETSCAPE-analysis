"""
  macro for analysis of jetscape events
  """

# This script allows one to run and analyze Jetscape events for a set of pT-hat bins.
# The user should pass in a yaml configuration file to specify the pT-hat binning and output directory.
#
# There are three main options, which should be run sequentially:
#     --run         Generate jetscape events by calling runJetscape
#     --analyze     Analyze Jetscape output file into a ROOT file
#     --plot        Make plots from the ROOT file
# For example: python doJetscapeAnalysis.py -c myConfig.yaml --run
#
# Requirements: (automatically satisfied by the Jetscape docker container)
#     * python 3.5 or higher
#     * ROOT (with pyroot support)
#     * If using --run option, requires Jetscape environment
#     * If using --analyze option, requires HepMC
#
# Author: James Mulligan (james.mulligan@berkeley.edu)

# General
import os
import sys
import argparse
import yaml
import subprocess
import fileinput
import shutil

# Data analysis and plotting
import ROOT
import common_base
import scaleHistograms
import plotResults

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
class doJetscapeAnalysis(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', xml_user_file='', xml_master_file='', run = False, analyze = False, plot = False, **kwargs):
    
    super(doJetscapeAnalysis, self).__init__(**kwargs)
    self.config_file = config_file
    self.xml_user_file = xml_user_file
    self.xml_master_file = xml_master_file
  
    self.run = run
    self.analyze = analyze
    self.plot = plot
    
    self.initialize_config()
    
    print(self)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
  
    # Read config file and create output directory, if needed
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Create output dir
    self.output_dir = config['output_dir']
    if not self.output_dir.endswith('/'):
      self.output_dir = self.output_dir + '/'
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

    # Get pT-hat bin list
    self.pt_hat_bins = config['pt_hat_bins']

    if self.plot:
      self.output_dir_pp = config['output_dir_pp']
      self.output_dir_AA = config['output_dir_AA']
      self.nEvents_pp = config['nEvents_pp']
      self.nEvents_AA = config['nEvents_AA']
      self.file_format = config['file_format']

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def do_jetscape_analysis(self):

    if self.run:   # Must be run inside Jetscape environment
      print("Run Jetscape events for all pT-hat bins!")
      self.runJetscape()

    if self.analyze:   # Must be run in environment with: HepMC, ROOT
      print("Analyze Jetscape output for all pT-hat bins!")
      self.analyzeOutput()
    
    if self.plot:   # Must be run in environment with: ROOT
      print("Plot histograms!")
      self.plotAnalysis()

  # ---------------------------------------------------------
  def runJetscape(self):

    # Loop through pT-hat bins and create directory structure with XML files for each bin
    for bin, pt_hat_min in enumerate(self.pt_hat_bins):
      
      # Set min,max of pT-hat bin
      if bin < ( len(self.pt_hat_bins) -1 ):
        pt_hat_max = self.pt_hat_bins[bin+1]
      else:
        continue
      
      # Create outputDir for each bin
      output_dir_bin = '{}{}'.format(self.output_dir, bin)
      if not output_dir_bin.endswith('/'):
        output_dir_bin = output_dir_bin + '/'
      if not os.path.exists(output_dir_bin):
        os.makedirs(output_dir_bin)

      # Set pT-hat values in Jetscape User XML configuration
      for line in fileinput.input(self.xml_user_file, inplace=True):
        if 'pTHatMin' in line:
          print('      <pTHatMin>{}</pTHatMin>'.format(pt_hat_min))
        elif 'pTHatMax' in line:
          print('      <pTHatMax>{}</pTHatMax>'.format(pt_hat_max))
        else:
          print(line, end='')
      shutil.copyfile(self.xml_user_file, '{}{}'.format(output_dir_bin, 'jetscape_user.xml'))
      shutil.copyfile(self.xml_master_file, '{}{}'.format(output_dir_bin, 'jetscape_master.xml'))

    # Loop through pt-hat bins and call Jetscape executable, and write output to pT-hat bin directory
    for bin, pt_hat_min in enumerate(self.pt_hat_bins):

      # Only iterate over lower bin edges
      if bin < ( len(self.pt_hat_bins) -1 ):
        pt_hat_max = self.pt_hat_bins[bin+1]
        print('PtHat: {} - {}'.format(pt_hat_min, pt_hat_max))
      else:
        continue

      # cd into bin directory in order to write Jetscape output there
      output_dir_bin = '{}{}'.format(self.output_dir, bin)
      os.chdir(output_dir_bin)

      # Run Jetscape executable
      logfileName = os.path.join(output_dir_bin, 'log_{}.txt'.format(bin))
      with open(logfileName, 'w') as logfile:
        cmd = '/home/jetscape-user/JETSCAPE/build/runJetscape jetscape_user.xml jetscape_master.xml'
        subprocess.run(cmd, check=True, shell=True, stdout=logfile)
      os.chdir(self.output_dir)

  # ---------------------------------------------------------
  def analyzeOutput(self):
    
    # Loop through pT-hat bins
    for bin, pt_hat_min in enumerate(self.pt_hat_bins):
      
      # Set min,max of pT-hat bin
      if bin < ( len(self.pt_hat_bins) -1 ):
        pt_hat_max = self.pt_hat_bins[bin+1]
        print('PtHat: {} - {}'.format(pt_hat_min, pt_hat_max))
      else:
        continue
      
      # Get outputDir for each bin
      output_dir_bin = '{}{}'.format(self.output_dir, bin)
      if not output_dir_bin.endswith('/'):
        output_dir_bin = output_dir_bin + '/'
      if not os.path.exists(output_dir_bin):
        print('output_dir_bin {} does not exist!'.format(bin))
    
      # Read HepMC output, get hadrons, do jet finding, and write histograms to ROOT file
      cmd = '/home/jetscape-user/JETSCAPE-analysis/build/runJetscapeAnalysis {} {}'.format(bin, output_dir_bin)
      subprocess.run(cmd, check=True, shell=True)

      # Scale histograms according to pthard bins cross-section
      print('Scaling pt-hat bins...')
      scaleHistograms.scaleHistograms(output_dir_bin, bin)
    
    # Merge all pthard bins into a single output file
    cmd = 'hadd {}AnalysisResultsFinal.root {}*/AnalysisResults.root'.format(self.output_dir, self.output_dir)
    subprocess.run(cmd, check=True, shell=True)

  # ---------------------------------------------------------
  def plotAnalysis(self):

    print('Plotting histograms...')
    plotResults.plotResults(self.output_dir_pp, self.output_dir_AA, self.nEvents_pp, self.nEvents_AA, self.file_format)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description="Plot analysis histograms")
  parser.add_argument("-c", "--configFile", action="store",
                      type=str, metavar="configFile",
                      default="config.yaml",
                      help="Path of config file for jetscape analysis")
  parser.add_argument("-u", "--xmlUserFile", action="store",
                      type=str, metavar="xmlUserFile",
                      default="/home/jetscape-user/JETSCAPE/config/jetscape_user.xml",
                      help="Path of JETSCAPE XML user file for jetscape analysis")
  parser.add_argument("-m", "--xmlMasterFile", action="store",
                      type=str, metavar="xmlMasterFile",
                      default="/home/jetscape-user/JETSCAPE/config/jetscape_master.xml",
                      help="Path of JETSCAPE XML master file for jetscape analysis")
  parser.add_argument("--run", action="store_true",
                      help="Whether to launch running of jetscape events")
  parser.add_argument("--analyze", action="store_true",
                      help="Whether to analyze output of jetscape events")
  parser.add_argument("--plot", action="store_true",
                      help="Whether to plot output of jetscape analysis")

  # Parse the arguments
  args = parser.parse_args()  
  print("Configuring...")

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print("File \"{0}\" does not exist! Exiting!".format(args.configFile))
    sys.exit(0)

  if args.run:
    print("XML user file: \"{0}\"".format(args.xmlUserFile))
    print("XML master file: \"{0}\"".format(args.xmlMasterFile))
    if not os.path.exists(args.xmlUserFile):
      print("File \"{0}\" does not exist! Exiting!".format(args.xmlUserFile))
      sys.exit(0)
    if not os.path.exists(args.xmlMasterFile):
      print("File \"{0}\" does not exist! Exiting!".format(args.xmlMasterFile))
      sys.exit(0)

  analysis = doJetscapeAnalysis(config_file = args.configFile, xml_user_file = args.xmlUserFile, xml_master_file = args.xmlMasterFile, run = args.run, analyze = args.analyze, plot = args.plot)
  analysis.do_jetscape_analysis()
