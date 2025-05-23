"""
  utils for histogramming / plotting

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import sys
import yaml
import argparse

# Data analysis and plotting
import ROOT
import ctypes
import numpy as np

# Base class
sys.path.append('.')
from jetscape_analysis.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotUtils(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # ---------------------------------------------------------------
    # Get bin array specified in config block
    # ---------------------------------------------------------------
    def bins_from_config(self, block, sqrts, observable_type, observable, centrality, centrality_index, suffix=''):

        if 'bins' in block:
            print(f'  Histogram with custom binning found for {observable} {centrality} {suffix}')
            return np.array(block['bins'])
        elif f'hepdata' in block:
            print(f'  Histogram with hepdata binning found for {observable} {centrality} {suffix}')
            return self.bins_from_hepdata(block, sqrts, observable_type, observable, centrality_index, suffix)
        else:
            print(f'  Warning: No binning found for {observable} {centrality} {suffix}')
            return np.array([])

    # ---------------------------------------------------------------
    # Get bin array from hepdata file specified in config block
    # ---------------------------------------------------------------
    def bins_from_hepdata(self, block, sqrts, observable_type, observable, centrality_index, suffix=''):

        # Open the HEPData file
        hepdata_dir = f'data/STAT/{sqrts}/{observable_type}/{observable}'
        hepdata_filename = os.path.join(hepdata_dir, block['hepdata'])
        f = ROOT.TFile(hepdata_filename, 'READ')

        # Find the relevant directory:
        # - We require that the centrality index is specified -- but
        #   the directory name may be found in the config either
        #   as a list of dirs, or a single dir with a list of histogram names
        # - The list of dir/hist names may also contain a suffix,
        #   which specifies e.g. the pt bin, jetR, or other parameters

        # First, check for dir names in config
        if f'hepdata_AA_dir{suffix}' in block:
            dir_key = f'hepdata_AA_dir{suffix}'
        elif f'hepdata_AA_dir' in block:
            dir_key = f'hepdata_AA_dir'
        else:
            print(f'hepdata_AA_dir{suffix} not found!')

        # Check for hist names in config
        if f'hepdata_AA_hname{suffix}' in block:
            h_key = f'hepdata_AA_hname{suffix}'
        elif f'hepdata_AA_hname' in block:
            h_key = f'hepdata_AA_hname'
        else:
            print(f'hepdata_AA_hname{suffix} not found!')
            return np.array([])

        # Get the appropriate centrality entry in the dir/hist list
        if type(block[dir_key]) is list:

            # If fewer entries than the observable's centrality, skip
            if centrality_index > len(block[dir_key])-1:
                return np.array([])

            dir_name = block[dir_key][centrality_index]
            h_name = block[h_key]

        else:

            # If fewer entries than the observable's centrality, skip
            if centrality_index > len(block[h_key])-1:
                return np.array([])

            dir_name = block[dir_key]
            h_name = block[h_key][centrality_index]

        # Get the histogram, and return the bins
        dir = f.Get(dir_name)
        h = dir.Get(h_name)
        bins = np.array(h.GetXaxis().GetXbins())

        # For certain Soft Drop observables, we need to exclude the "untagged" bin so that it will become underflow
        if observable in ["zg_alice", "tg_alice"]:
            bins = bins[1:]
        # For ATLAS RAA y-dependence, add a bin at 0<abs(y)<0.3 (it is not included in the hepdata since they take a ratio to that bin)
        if observable == 'pt_y_atlas':
            bins = np.insert(bins, 0, 0., axis=0)

        f.Close()

        return bins

    # ---------------------------------------------------------------
    # Get tgraph from hepdata file specified in config block
    # ---------------------------------------------------------------
    def tgraph_from_hepdata(self, block, is_AA, sqrts, observable_type, observable, centrality_index, suffix='', pt_suffix=''):

        # Open the HEPData file
        hepdata_dir = f'data/STAT/{sqrts}/{observable_type}/{observable}'
        hepdata_filename = os.path.join(hepdata_dir, block['hepdata'])
        f = ROOT.TFile(hepdata_filename, 'READ')

        # Find the relevant directory:
        # - The list of dir/hist names may contain a suffix,
        #   which specifies e.g. the pt bin, jetR, or other parameters

        if is_AA:
            system = 'AA'
        else:
            system = 'pp'

        # First, check for dir names in config
        if f'hepdata_{system}_dir{suffix}' in block:
            dir_key = f'hepdata_{system}_dir{suffix}'
        elif f'hepdata_{system}_dir{suffix}{pt_suffix}' in block:
            dir_key = f'hepdata_{system}_dir{suffix}{pt_suffix}'
        elif f'hepdata_{system}_dir' in block:
            dir_key = f'hepdata_{system}_dir'
        else:
            #print(f'hepdata_pp_dir{suffix} not found!')
            return None

        # Check for hist names in config
        if f'hepdata_{system}_gname{suffix}' in block:
            g_key = f'hepdata_{system}_gname{suffix}'
        elif f'hepdata_{system}_gname{suffix}{pt_suffix}' in block:
            g_key = f'hepdata_{system}_gname{suffix}{pt_suffix}'
        elif f'hepdata_{system}_gname' in block:
            g_key = f'hepdata_{system}_gname'
        else:
            #print(f'hepdata_{system}_gname{suffix} not found!')
            return None

        # Get the appropriate centrality entry in the dir/hist list
        if type(block[dir_key]) is list:

            # If fewer entries than the observable's centrality, skip
            if centrality_index > len(block[dir_key])-1:
                return np.array([])

            dir_name = block[dir_key][centrality_index]
            g_name = block[g_key]

        elif type(block[g_key]) is list:

            # If fewer entries than the observable's centrality, skip
            if centrality_index > len(block[g_key])-1:
                return np.array([])

            dir_name = block[dir_key]
            g_name = block[g_key][centrality_index]

        else:

            dir_name = block[dir_key]
            g_name = block[g_key]

        # Get the tgraph, and return the bins
        dir = f.Get(dir_name)
        g = dir.Get(g_name)
        f.Close()

        return g

    # ---------------------------------------------------------------
    # Get tgraph from data points specified in config block
    # ---------------------------------------------------------------
    def tgraph_from_yaml(self, block, is_AA, sqrts, observable_type, observable, centrality_index, suffix='', pt_suffix=''):

        # Load yaml file containing the data
        data_dir = f'data/STAT/{sqrts}/{observable_type}/{observable}'
        data_filename = os.path.join(data_dir, block['custom_data'])
        with open(data_filename, 'r') as stream:
            data = yaml.safe_load(stream)

        # Find the relevant config block:
        # - The list of dir/hist names may contain a suffix,
        #   which specifies e.g. the pt bin, jetR, or other parameters

        if is_AA:
            system = 'AA'
        else:
            system = 'pp'

        # First, check for key names in config
        if f'data_{system}{suffix}' in data:
            key = f'data_{system}{suffix}'
        elif f'data_{system}{suffix}{pt_suffix}' in data:
            key = f'data_{system}{suffix}{pt_suffix}'
        elif f'data_{system}' in data:
            key = f'data_{system}'
        else:
            #print(f'data_pp_dir{suffix} not found!')
            return None

        # Get the appropriate centrality entry in the list
        if 'x' in data[key]:
            x = np.array(data[key]['x'][centrality_index])
        else:
            bins = np.array(block['bins'])
            x = (bins[1:] + bins[:-1]) / 2
        n = len(x)

        if isinstance(data[key]['y'][0], list):
            if len(data[key]['y']) > centrality_index:
                y = np.array(data[key]['y'][centrality_index])
                if 'y_err' in data[key]:
                    y_err = np.array(data[key]['y_err'][centrality_index])
                else:
                    y_err = np.zeros(n)
            else:
                return None
        else:
            y = np.array(data[key]['y'])
            if 'y_err' in data[key]:
                y_err = np.array(data[key]['y_err'])
            else:
                y_err = np.zeros(n)

        # Construct TGraph
        g = ROOT.TGraphAsymmErrors(n, x, y, np.zeros(n), np.zeros(n), y_err, y_err)

        return g

    #---------------------------------------------------------------
    # Truncate data tgraph to histogram binning range
    #---------------------------------------------------------------
    def truncate_tgraph(self, g, h, is_AA=False):

        #print('truncate_tgraph')
        #print(h.GetName())
        #print(np.array(h.GetXaxis().GetXbins()))

        # Create new TGraph with number of points equal to number of histogram bins
        nBins = h.GetNbinsX()
        g_new = ROOT.TGraphAsymmErrors(nBins)
        g_new.SetName(f'{g.GetName()}_truncated')

        h_offset = 0
        for bin in range(1, nBins+1):

            # Get histogram (x)
            h_x = h.GetBinCenter(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'gy: {gy}')

            # If traph is offset from center of the bin, center it
            xErrLow = g.GetErrorXlow(bin-1)
            xErrUp = g.GetErrorXhigh(bin-1)
            if xErrLow > 0 and xErrUp > 0:
                x_min = gx - xErrLow
                x_max = gx + xErrUp
                x_center = (x_min + x_max)/2.
                if h_x > x_min and h_x < x_max:
                    if not np.isclose(gx, x_center):
                        gx = x_center

            # If tgraph starts below hist (e.g. when hist has min cut), try to get next tgraph point
            g_offset = 0
            while gx+1e-8 < h_x and g_offset < g.GetN()+1:
                g_offset += 1
                gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1+g_offset)
            #print(f'new gx: {gx}')

            # If tgraph started above hist (see below) and we exhausted the tgraph points, skip
            if h_offset > 0 and np.isclose(gx, 0):
                continue

            # If tgraph starts above hist, try to get next hist bin
            h_offset = 0
            while gx-1e-8 > h_x and h_offset < nBins+1:
                h_offset += 1
                h_x = h.GetBinCenter(bin+h_offset)
                #print(f'h_x: {h_x}')
                #print(f'gx: {gx}')

            if not np.isclose(h_x, gx):
                if is_AA:
                    sys.exit(f'ERROR: hist x: {h_x}, graph x: {gx}')
                else:
                    print(f'WARNING: hist x: {h_x}, graph x: {gx}')
                    return None

            g_new.SetPoint(bin-1, gx, gy)
            g_new.SetPointError(bin-1, 0, 0, yErrLow, yErrUp)
            #print()

        return g_new

    #---------------------------------------------------------------
    # Divide a histogram by a tgraph, point-by-point
    #---------------------------------------------------------------
    def divide_histogram_by_tgraph(self, h, g, include_tgraph_uncertainties=True):

        # Truncate tgraph to range of histogram bins
        g_truncated = self.truncate_tgraph(g, h)
        if not g_truncated:
            return None

        # Clone tgraph, in order to return a new one
        g_new = g_truncated.Clone(f'{g_truncated.GetName()}_divided')

        nBins = h.GetNbinsX()
        for bin in range(1, nBins+1):

            # Get histogram (x,y)
            h_x = h.GetBinCenter(bin)
            h_y = h.GetBinContent(bin)
            h_error = h.GetBinError(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g_truncated, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'h_y: {h_y}')
            #print(f'gy: {gy}')

            if not np.isclose(h_x, gx):
                print(f'WARNING: hist x: {h_x}, graph x: {gx} -- will not plot ratio')
                return None

            new_content = h_y / gy

            # Combine tgraph and histogram relative uncertainties in quadrature
            if gy > 0. and h_y > 0.:
                if include_tgraph_uncertainties:
                    new_error_low = np.sqrt( pow(yErrLow/gy,2) + pow(h_error/h_y,2) ) * new_content
                    new_error_up = np.sqrt( pow(yErrUp/gy,2) + pow(h_error/h_y,2) ) * new_content
                else:
                    new_error_low = h_error/h_y * new_content
                    new_error_up = h_error/h_y * new_content
            else:
                new_error_low = (yErrLow/gy) * new_content
                new_error_up = (yErrUp/gy) * new_content

            g_new.SetPoint(bin-1, h_x, new_content)
            g_new.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)

        return g_new

    #---------------------------------------------------------------
    # Get points from tgraph by index
    #---------------------------------------------------------------
    def get_gx_gy(self, g, index):

        g_x = ctypes.c_double(0)
        g_y = ctypes.c_double(0)
        g.GetPoint(index, g_x, g_y)
        yErrLow = g.GetErrorYlow(index)
        yErrUp  = g.GetErrorYhigh(index)

        gx = g_x.value
        gy = g_y.value

        return gx, gy, yErrLow, yErrUp

    #---------------------------------------------------------------
    # Divide a tgraph by a tgraph, point-by-point: g1/g2
    # NOTE: Ignore uncertainties on denominator
    #---------------------------------------------------------------
    def divide_tgraph_by_tgraph(self, g1, g2):

        # Clone tgraph, in order to return a new one
        g_new = g1.Clone(f'{g1.GetName()}_divided')

        if g1.GetN() != g2.GetN():
            sys.exit(f'ERROR: TGraph {g1.GetName()} has {g1.GetN()} points, but {g2.GetName()} has {g2.GetN()} points')

        for i in range(0, g1.GetN()):

            # Get TGraph (x,y) and errors
            g1_x = ctypes.c_double(0)
            g1_y = ctypes.c_double(0)
            g1.GetPoint(i, g1_x, g1_y)
            y1ErrLow = g1.GetErrorYlow(i)
            y1ErrUp  = g1.GetErrorYhigh(i)
            g1x = g1_x.value
            g1y = g1_y.value

            g2_x = ctypes.c_double(0)
            g2_y = ctypes.c_double(0)
            g2.GetPoint(i, g2_x, g2_y)
            g2x = g2_x.value
            g2y = g2_y.value

            if not np.isclose(g1x, g2x):
                sys.exit(f'ERROR: TGraph {g1.GetName()} point {i} at {g1x}, but {g2.GetName()} at {g2x}')

            new_content = g1y / g2y
            new_error_low = y1ErrLow/g1y * new_content
            new_error_up = y1ErrUp/g1y * new_content

            g_new.SetPoint(i, g1x, new_content)
            g_new.SetPointError(i, 0, 0, new_error_low, new_error_up)
        return g_new

    #-------------------------------------------------------------------------------------------
    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def setup_legend(self, leg, textSize, sep, title=''):

        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        leg.SetMargin(0.25)
        leg.SetTextSize(textSize)
        leg.SetEntrySeparation(sep)
        leg.SetHeader(title)

    #-------------------------------------------------------------------------------------------
    def setOptions(self):

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
