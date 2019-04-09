"""
  macro for analysis of jetscape events in hepmc format
  """

# This script allows one to run and analyze Jetscape events for a set of pT-hat bins.
# The user should pass in a yaml configuration file to specify the pT-hat binning and output directory.
#
# There are three main options, which should be run sequentially:
#     --run         Generate jetscape events by calling runJetscape
#     --analyze     Analyze Jetscape output file into a ROOT file
#     --plot        Make plots from the ROOT file
# For example: python doJetAnalysis.py -c myConfig.yaml --run
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

# Data analysis and plotting
import ROOT

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def doJetAnalysis(configFile, run, analyze, plot, fileFormat):
  
  # Read config file and create output directory, if needed
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)
  
  # Create output dir
  outputDir = config['outputDir']
  if not outputDir.endswith('/'):
    outputDir = outputDir + '/'
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  # Get pT-hat bin list
  PtHatBins = config['PtHatBins']
  print('Pt-hat bins: {}'.format(PtHatBins))

  #############################################

  if run:   # Must be run inside Jetscape environment
    print("Run Jetscape events for all pT-hat bins!")
    runJetscape(PtHatBins, outputDir, fileFormat)

  if analyze:   # Must be run in environment with: HepMC, ROOT
    print("Analyze Jetscape output for all pT-hat bins!")
    analyzeOutput(PtHatBins, outputDir, fileFormat)
    
  if plot:   # Must be run in environment with: ROOT
    print("Plot histograms!")
    plotAnalysis(outputDir, fileFormat)

# ---------------------------------------------------------
def runJetscape(PtHatBins, outputDir, fileFormat):

  # Loop through pT-hat bins
  for bin, PtHatMin in enumerate(PtHatBins):
    
    # Set min,max of pT-hat bin
    if bin < ( len(PtHatBins) -1 ):
      PtHatMax = PtHatBins[bin+1]
      print('PtHat: {} - {}'.format(PtHatMin, PtHatMax))
    else:
      continue
    
    # Create outputDir for each bin
    outputDirBin = '{}{}'.format(outputDir, bin)
    if not outputDirBin.endswith('/'):
      outputDirBin = outputDirBin + '/'
    if not os.path.exists(outputDirBin):
      os.makedirs(outputDirBin)

    # Call Jetscape executable
    logfileName = os.path.join(outputDirBin, 'log_{}.txt'.format(bin))
    with open(logfileName, 'w') as logfile:
      cmd = 'cd /home/jetscape-user/JETSCAPE/build && ./runJetscape'
      subprocess.run(cmd, check=True, shell=True, stdout=logfile)

    # Move output file into pT-hat bin directory
    cmd = 'mv /home/jetscape-user/JETSCAPE/build/test_out.hepmc {}'.format(outputDirBin)
    subprocess.run(cmd, check=True, shell=True)

# ---------------------------------------------------------
def analyzeOutput(PtHatBins, outputDir, fileFormat):
  
  # Loop through pT-hat bins
  for bin, PtHatMin in enumerate(PtHatBins):
    
    # Set min,max of pT-hat bin
    if bin < ( len(PtHatBins) -1 ):
      PtHatMax = PtHatBins[bin+1]
      print('PtHat: {} - {}'.format(PtHatMin, PtHatMax))
    else:
      continue
  
    cmd = '../JETSCAPE/build/readerTestHepMC test_out.hepmc'
    #subprocess.run(cmd, check=True, shell=True)

# ---------------------------------------------------------
def plotAnalysis(outputDir, fileFormat):

  print('plot')

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
  parser.add_argument("--run", action="store_true",
                      help="Whether to launch running of jetscape events")
  parser.add_argument("--analyze", action="store_true",
                      help="Whether to analyze output of jetscape events")
  parser.add_argument("--plot", action="store_true",
                      help="Whether to plot output of jetscape analysis")
  parser.add_argument("-i", "--imageFormat", action="store",
                      type=str, metavar="imageFormat",
                      default=".pdf",
                      help="Image format to save plots in, e.g. \".pdf\" or \".png\"")
  
  # Parse the arguments
  args = parser.parse_args()
  
  print("Configuring...")
  print("configFile: \"{0}\"".format(args.configFile))
  print("imageFormat: \"{0}\"".format(args.imageFormat))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print("File \"{0}\" does not exist! Exiting!".format(args.configFile))
    sys.exit(0)

doJetAnalysis(configFile = args.configFile, run = args.run, analyze = args.analyze, plot = args.plot, fileFormat = args.imageFormat)
