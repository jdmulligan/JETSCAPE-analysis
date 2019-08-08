"""
  macro for plotting analyzed jetscape events
  """

# This script plots histograms created in the analysis of Jetscape events
#
# Author: James Mulligan (james.mulligan@berkeley.edu)

# General
import os
import sys

# Data analysis and plotting
import ROOT
from array import *

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
def plotResults(outputDir_pp, outputDir_AA, nEvents_pp, nEvents_AA, fileFormat):
  
  setOptions()
  ROOT.gROOT.ForceStyle()
  
  #plotPPJetCrossSection(outputDir_pp, nEvents_pp, fileFormat)
  plotPPJetCrossSectionAmit(outputDir_pp, nEvents_pp, fileFormat)
  #plotHadronRAA(outputDir_pp, outputDir_AA, nEvents_pp, nEvents_AA, fileFormat)

#-------------------------------------------------------------------------------------------
def plotPPJetCrossSection(outputDir_pp, nEvents_pp, fileFormat):

  filename_pp = '{}/AnalysisResultsFinal.root'.format(outputDir_pp)
  f_pp = ROOT.TFile(filename_pp, 'READ')
  
  hJetPt_pp = f_pp.Get('hJet040_PtScaled')
  hJetPt_pp.Rebin(100)
  hJetPt_pp.Scale(1/nEvents_pp, 'width')
  
  #radius = label[1:-5]      ##cut the string into place to find the right radius
  #radius = 0.1*int(radius)  ##turn the string number into a integer
  #etaAcc = 0.7 - radius

  # Create canvas
  cSpectrum = ROOT.TCanvas("cSpectrum","cSpectrum: hist",600,600)
  cSpectrum.Draw()
  cSpectrum.cd()
  cSpectrum.SetLogy()

  # Set pad and histo arrangement
  myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
  myPad.SetLeftMargin(0.22)
  myPad.SetTopMargin(0.04)
  myPad.SetRightMargin(0.04)
  myPad.SetBottomMargin(0.15)
  
  myBlankHisto = ROOT.TH1F("myBlankHisto","Blank Histogram",20,0,200)
  myBlankHisto.SetNdivisions(505)
  myBlankHisto.SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
  myBlankHisto.GetYaxis().SetTitleOffset(2.2)
  myBlankHisto.SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T,jet}d#it{#eta}_{jet}} #left[mb (GeV/c)^{-1}#right]")
  myBlankHisto.SetMaximum(9e-3)
  myBlankHisto.SetMinimum(2e-9)
  myBlankHisto.Draw()

  # Set spectra styles
  print("Plot pp spectrum")
  hJetPt_pp.SetMarkerSize(2)
  hJetPt_pp.SetMarkerStyle(33)
  hJetPt_pp.SetMarkerColor(600-6)
  hJetPt_pp.SetLineStyle(1)
  hJetPt_pp.SetLineWidth(2)
  hJetPt_pp.SetLineColor(600-6)
  
  # Draw spectra
  hJetPt_pp.Draw('PE X0 same')

  # # # # # # # # # # # # # # # # # # # # # # # #
  # Add legends and text
  # # # # # # # # # # # # # # # # # # # # # # # #
  system = ROOT.TLatex(0.49,0.835,"pp  #sqrt{#it{s}} = 2.76 TeV")
  system.SetNDC()
  system.SetTextSize(0.044)
  system.Draw()
  system2 = ROOT.TLatex(0.49,0.90,"JETSCAPE")
  system2.SetNDC()
  system2.SetTextSize(0.044)
  system2.Draw()
  
  myLegend = ROOT.TLegend(0.4,0.7,0.7,0.82)
  setupLegend(myLegend, 0.035*1.1) #0.04
  myLegend.AddEntry(myBlankHisto,"Anti-#it{k}_{T} #it{R} = 0.2 | #it{#eta}_{jet}| < 1.0","")
  myLegend.Draw()
  
  myLegend2pp = ROOT.TLegend(0.65,0.6,0.8,0.68)
  setupLegend(myLegend2pp,0.04)
  myLegend2pp.AddEntry(hJetPt_pp,"pp","Pe")
  myLegend2pp.Draw()
  
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Save
  # # # # # # # # # # # # # # # # # # # # # # # #
  outputFilename = os.path.join(outputDir_pp, "hJetCrossSectionPP{}".format(fileFormat))
  cSpectrum.SaveAs(outputFilename)

#-------------------------------------------------------------------------------------------
def plotPPJetCrossSectionAmit(outputDir_pp, nEvents_pp, fileFormat):
  
  filename_pp = '{}/AnalysisResultsFinal.root'.format(outputDir_pp)
  f_pp = ROOT.TFile(filename_pp, 'READ')
  
  #amitPt = ([75, 85, 95, 105, 120, 140, 160, 180, 200, 225, 255, 285])
  #binArray = ([70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300])
  binArray = ([70, 80, 90, 100, 110, 130, 150, 170, 190])
  amitCrossSection = [3.73703e-06, 1.82665e-06, 9.5494e-07, 5.3191e-07, 2.44257e-07, 9.25817e-08, 3.86036e-08, 1.75827e-08, 8.47066e-09, 3.70068e-09, 1.44942e-09, 6.09325e-10]
  amitUncertainty = [2.5731e-08, 8.03882e-09, 4.37744e-09, 2.94087e-09, 8.87617e-10, 3.34769e-10, 1.39597e-10, 6.11302e-11, 2.98531e-11, 1.04712e-11, 4.38574e-12, 1.82346e-12]
  
  nBins = len(binArray) - 1
  bin_array = array('d',binArray)
  
  hJetPt_pp_original = f_pp.Get('hJet040_PtScaled')
  #hJetPt_pp.Rebin(100)
  hJetPt_pp = hJetPt_pp_original.Rebin(nBins, '{}_NewBinning'.format(hJetPt_pp_original.GetName()), bin_array)
  hJetPt_pp.Scale(1/nEvents_pp, 'width')
  
  hJetPt_pp_amit = hJetPt_pp.Clone()
  hJetPt_pp_amit.SetName('hJetPt_pp_amit')
  for bin in range(len(binArray)-1):
    print('{} : {}: {}'.format(bin+1, hJetPt_pp_amit.GetBinCenter(bin+1), amitCrossSection[bin]))
    hJetPt_pp_amit.SetBinContent(bin+1, amitCrossSection[bin])
  hJetPt_pp_amit.SetMarkerSize(1)
  hJetPt_pp_amit.SetMarkerStyle(21)
  hJetPt_pp_amit.SetMarkerColor(1)

  #radius = label[1:-5]      ##cut the string into place to find the right radius
  #radius = 0.1*int(radius)  ##turn the string number into a integer
  #etaAcc = 0.7 - radius

  # Create canvas
  c = ROOT.TCanvas("c","c: pT",800,850)
  c.cd()
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0)
  pad1.SetLeftMargin(0.15)
  pad1.SetRightMargin(0.05)
  pad1.SetTopMargin(0.05)
  pad1.SetLogy()
  pad1.Draw()
  pad1.cd()
  
  # Set pad and histo arrangement
  myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
  myPad.SetLeftMargin(0.22)
  myPad.SetTopMargin(0.04)
  myPad.SetRightMargin(0.04)
  myPad.SetBottomMargin(0.15)
  
  myBlankHisto = ROOT.TH1F("myBlankHisto","Blank Histogram", 120, 70, 190)
  myBlankHisto = myBlankHisto.Rebin(nBins, '{}_NewBinning'.format(myBlankHisto.GetName()), bin_array)

  myBlankHisto.SetNdivisions(505)
  myBlankHisto.SetXTitle("#it{p}_{T,jet} (GeV/#it{c})")
  myBlankHisto.GetYaxis().SetTitleOffset(2.2)
  myBlankHisto.SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T,jet}d#it{#eta}_{jet}} #left[mb (GeV/c)^{-1}#right]")
  myBlankHisto.SetMaximum(9e-3)
  myBlankHisto.SetMinimum(2e-9)
  myBlankHisto.Draw()
  
  # Set spectra styles
  print("Plot pp spectrum")
  hJetPt_pp.SetMarkerSize(2)
  hJetPt_pp.SetMarkerStyle(33)
  hJetPt_pp.SetMarkerColor(600-6)
  hJetPt_pp.SetLineStyle(1)
  hJetPt_pp.SetLineWidth(2)
  hJetPt_pp.SetLineColor(600-6)
  
  # Draw spectra
  hJetPt_pp.Draw('PE X0 same')
  hJetPt_pp_amit.Draw('PE X0 same')

  # # # # # # # # # # # # # # # # # # # # # # # #
  # Add legends and text
  # # # # # # # # # # # # # # # # # # # # # # # #
  system = ROOT.TLatex(0.49,0.835,"pp  #sqrt{#it{s}} = 2.76 TeV")
  system.SetNDC()
  system.SetTextSize(0.044)
  system.Draw()
  system2 = ROOT.TLatex(0.49,0.90,"JETSCAPE")
  system2.SetNDC()
  system2.SetTextSize(0.044)
  system2.Draw()
  
  myLegend = ROOT.TLegend(0.4,0.7,0.7,0.82)
  setupLegend(myLegend, 0.035*1.1) #0.04
  myLegend.AddEntry(myBlankHisto,"Anti-#it{k}_{T} #it{R} = 0.2 | #it{#eta}_{jet}| < 1.0","")
  myLegend.Draw()
  
  myLegend2pp = ROOT.TLegend(0.65,0.6,0.8,0.68)
  setupLegend(myLegend2pp,0.04)
  myLegend2pp.AddEntry(hJetPt_pp,"Me","Pe")
  myLegend2pp.AddEntry(hJetPt_pp_amit,"Amit","Pe")
  myLegend2pp.Draw()

  c.cd()
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.35)
  pad2.SetLeftMargin(0.15)
  pad2.SetRightMargin(0.05)
  pad2.Draw()
  pad2.cd()

  # plot ratio
  hRatio = hJetPt_pp.Clone()
  hRatio.Divide(hJetPt_pp_amit)
  hRatio.SetMarkerStyle(21)
  hRatio.SetMarkerColor(1)

  hRatio.GetXaxis().SetRangeUser(binArray[0], binArray[-1])
  hRatio.GetXaxis().SetTitleSize(30)
  hRatio.GetXaxis().SetTitleFont(43)
  hRatio.GetXaxis().SetTitleOffset(4.)
  hRatio.GetXaxis().SetLabelFont(43)
  hRatio.GetXaxis().SetLabelSize(20)
  hRatio.GetXaxis().SetTitle('pT')
  
  hRatio.GetYaxis().SetTitle('Me/Amit')
  hRatio.GetYaxis().SetTitleSize(20)
  hRatio.GetYaxis().SetTitleFont(43)
  hRatio.GetYaxis().SetTitleOffset(2.2)
  hRatio.GetYaxis().SetLabelFont(43)
  hRatio.GetYaxis().SetLabelSize(20)
  hRatio.GetYaxis().SetNdivisions(505)

  hRatio.Draw("P E")

  outputFilename = os.path.join(outputDir_pp, "hJetCrossSectionPP_Ratio{}".format(fileFormat))
  c.SaveAs(outputFilename)

#-------------------------------------------------------------------------------------------
def plotHadronRAA(outputDir_pp, outputDir_AA, nEvents_pp, nEvents_AA, fileFormat):

  filename_pp = '{}/AnalysisResultsFinal.root'.format(outputDir_pp)
  f_pp = ROOT.TFile(filename_pp, 'READ')
  hHadronPt_pp = f_pp.Get('hHadronPtScaled')
  hHadronPt_pp.SetName('hHadronPt_pp')
  hHadronPt_pp.Rebin(50)
  hHadronPt_pp.Scale(1/nEvents_pp, 'width')

  filename_AA = '{}/AnalysisResultsFinal.root'.format(outputDir_AA)
  f_AA = ROOT.TFile(filename_AA, 'READ')
  hHadronPt_AA = f_AA.Get('hHadronPtScaled')
  hHadronPt_AA.SetName('hHadronPt_AA')
  hHadronPt_AA.Rebin(50)
  hHadronPt_AA.Scale(1/nEvents_AA, 'width')

  hHadronRaa = hHadronPt_AA.Clone()
  hHadronRaa.SetName('hHadronRaa')
  hHadronRaa.Divide(hHadronPt_pp)

  # Create canvas
  cSpectrum = ROOT.TCanvas("cSpectrum","cSpectrum: hist",600,600)
  cSpectrum.Draw()
  cSpectrum.cd()

  # Set spectra styles
  print("Plot hadron RAA")
  hHadronRaa.SetMarkerSize(2)
  hHadronRaa.SetMarkerStyle(21)
  hHadronRaa.SetMarkerColor(600-6)
  
  # Draw spectra
  hHadronRaa.Draw("PE X0 same")

  outputFilename = os.path.join(outputDir_pp, "hHadronRAA{}".format(fileFormat))
  cSpectrum.SaveAs(outputFilename)

# Set legend parameters
#-------------------------------------------------------------------------------------------
def setupLegend(leg, textSize):
  
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetMargin(0.25);
  leg.SetTextSize(textSize);
  leg.SetEntrySeparation(0.5);

#-------------------------------------------------------------------------------------------
def setOptions():
  
  font = 42
  
  ROOT.gStyle.SetFrameBorderMode(0)
  ROOT.gStyle.SetFrameFillColor(0)
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetPadBorderMode(0)
  ROOT.gStyle.SetPadColor(10)
  ROOT.gStyle.SetCanvasColor(10)
  ROOT.gStyle.SetTitleFillColor(10)
  ROOT.gStyle.SetTitleBorderSize(1)
  ROOT.gStyle.SetStatColor(10)
  ROOT.gStyle.SetStatBorderSize(1)
  ROOT.gStyle.SetLegendBorderSize(1)
  
  ROOT.gStyle.SetDrawBorder(0)
  ROOT.gStyle.SetTextFont(font)
  ROOT.gStyle.SetStatFont(font)
  ROOT.gStyle.SetStatFontSize(0.05)
  ROOT.gStyle.SetStatX(0.97)
  ROOT.gStyle.SetStatY(0.98)
  ROOT.gStyle.SetStatH(0.03)
  ROOT.gStyle.SetStatW(0.3)
  ROOT.gStyle.SetTickLength(0.02,"y")
  ROOT.gStyle.SetEndErrorSize(3)
  ROOT.gStyle.SetLabelSize(0.05,"xyz")
  ROOT.gStyle.SetLabelFont(font,"xyz")
  ROOT.gStyle.SetLabelOffset(0.01,"xyz")
  ROOT.gStyle.SetTitleFont(font,"xyz")
  ROOT.gStyle.SetTitleOffset(1.2,"xyz")
  ROOT.gStyle.SetTitleSize(0.045,"xyz")
  ROOT.gStyle.SetMarkerSize(1)
  ROOT.gStyle.SetPalette(1)
  
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptFit(0)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
  print("Executing plotResults.py...")
  print("")
  
  plotResults(outputDir_pp, outputDir_AA, nEvents_pp, nEvents_AA, fileFormat)
