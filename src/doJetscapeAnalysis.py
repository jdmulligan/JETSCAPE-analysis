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
import fileinput
import shutil

# Data analysis and plotting
import ROOT
import scaleHistograms

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def doJetscapeAnalysis(configFile, xmlUserFile, xmlMasterFile, run, analyze, plot, fileFormat):
  
  # Read config file and create output directory, if needed
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)
  
  # Create output dir
  outputDir = config['outputDir']
  if not outputDir.endswith('/'):
    outputDir = outputDir + '/'
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)
  print('Output directory: {}'.format(outputDir))

  # Get pT-hat bin list
  PtHatBins = config['PtHatBins']
  print('Pt-hat bins: {}'.format(PtHatBins))

  #############################################

  if run:   # Must be run inside Jetscape environment
    print("Run Jetscape events for all pT-hat bins!")
    runJetscape(PtHatBins, xmlUserFile, xmlMasterFile, outputDir, fileFormat)

  if analyze:   # Must be run in environment with: HepMC, ROOT
    print("Analyze Jetscape output for all pT-hat bins!")
    analyzeOutput(PtHatBins, outputDir, fileFormat)
    
  if plot:   # Must be run in environment with: ROOT
    print("Plot histograms!")
    plotAnalysis(outputDir, fileFormat)

# ---------------------------------------------------------
def runJetscape(PtHatBins, xmlUserFile, xmlMasterFile, outputDir, fileFormat):

  # Loop through pT-hat bins and create directory structure with XML files for each bin
  for bin, PtHatMin in enumerate(PtHatBins):
    
    # Set min,max of pT-hat bin
    if bin < ( len(PtHatBins) -1 ):
      PtHatMax = PtHatBins[bin+1]
    else:
      continue
    
    # Create outputDir for each bin
    outputDirBin = '{}{}'.format(outputDir, bin)
    if not outputDirBin.endswith('/'):
      outputDirBin = outputDirBin + '/'
    if not os.path.exists(outputDirBin):
      os.makedirs(outputDirBin)

    # Set pT-hat values in Jetscape User XML configuration
    for line in fileinput.input(xmlUserFile, inplace=True):
      if 'pTHatMin' in line:
        print('      <pTHatMin>{}</pTHatMin>'.format(PtHatMin))
      elif 'pTHatMax' in line:
        print('      <pTHatMax>{}</pTHatMax>'.format(PtHatMax))
      else:
        print(line, end='')
    shutil.copyfile(xmlUserFile, '{}{}'.format(outputDirBin, 'jetscape_user.xml'))
    shutil.copyfile(xmlMasterFile, '{}{}'.format(outputDirBin, 'jetscape_master.xml'))

  # Loop through pt-hat bins and call Jetscape executable, and write output to pT-hat bin directory    
  for bin, PtHatMin in enumerate(PtHatBins):

    # Only iterate over lower bin edges
    if bin < ( len(PtHatBins) -1 ):
      PtHatMax = PtHatBins[bin+1]
      print('PtHat: {} - {}'.format(PtHatMin, PtHatMax))
    else:
      continue

    # cd into bin directory in order to write Jetscape output there
    outputDirBin = '{}{}'.format(outputDir, bin)
    os.chdir(outputDirBin)

    # Run Jetscape executable
    logfileName = os.path.join(outputDirBin, 'log_{}.txt'.format(bin))
    with open(logfileName, 'w') as logfile:
      cmd = '/home/jetscape-user/JETSCAPE/build/runJetscape {} {}'.format(xmlUserFile, xmlMasterFile)
      subprocess.run(cmd, check=True, shell=True, stdout=logfile)
    os.chdir(outputDir)

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
        
    # Get outputDir for each bin
    outputDirBin = '{}{}'.format(outputDir, bin)
    if not outputDirBin.endswith('/'):
      outputDirBin = outputDirBin + '/'
    if not os.path.exists(outputDirBin):
      print('outputDirBin {} does not exist!'.format(bin))
  
    # Read HepMC output, get hadrons, do jet finding, and write histograms to ROOT file
    cmd = '/home/jetscape-user/JETSCAPE-analysis/build/runJetscapeAnalysis {} {}'.format(bin, outputDirBin)
    subprocess.run(cmd, check=True, shell=True)

    # Scale histograms according to pthard bins cross-section
    print('Scaling pt-hat bins...')
    scaleHistograms.scaleHistograms(outputDirBin, bin)
            
  # Merge all pthard bins into a single output file
  cmd = 'hadd {}AnalysisResultsFinal.root {}*/AnalysisResults.root'.format(outputDir, outputDir)
  subprocess.run(cmd, check=True, shell=True)

# ---------------------------------------------------------
def plotAnalysis(outputDir, fileFormat):

  print('Plotting histograms...')

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
  parser.add_argument("-i", "--imageFormat", action="store",
                      type=str, metavar="imageFormat",
                      default=".pdf",
                      help="Image format to save plots in, e.g. \".pdf\" or \".png\"")
  
  # Parse the arguments
  args = parser.parse_args()  
  print("Configuring...")

  # If invalid configFile is given, exit
  print("configFile: \"{0}\"".format(args.configFile))
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

  if args.plot:
    print("imageFormat: \"{0}\"".format(args.imageFormat))
      
doJetscapeAnalysis(configFile = args.configFile, xmlUserFile = args.xmlUserFile, xmlMasterFile = args.xmlMasterFile, run = args.run, analyze = args.analyze, plot = args.plot, fileFormat = args.imageFormat)
