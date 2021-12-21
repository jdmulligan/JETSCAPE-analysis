# Quick script to get pp and AA distributions from TH1s/TGraphs and form PbPb/pp ratio

import ROOT
import numpy as np

#------------------------------------------
def tgraph_to_numpy(g):

    x = np.zeros(g.GetN())
    y = np.zeros(g.GetN())
    for i in range(g.GetN()):
        g_x = ROOT.Double(0)
        g_y = ROOT.Double(0)
        g.GetPoint(i, g_x, g_y)
        x[i] = g_x
        y[i] = g_y

    xerr = np.array([g.GetErrorX(i) for i in range(g.GetN())])
    yerr = np.array([g.GetErrorY(i) for i in range(g.GetN())])

    x_bins = np.append(x-xerr, x[-1]+xerr[-1])

    # Only include those points with nonzero y values
    truncate_index = np.argmax(y > 0.)

    return x_bins[truncate_index:], y[truncate_index:], yerr[truncate_index:]

#------------------------------------------
# Main

n_pt_bins = 4
n_cent_bins = 4

filename_pp = 'pp.R04.Aug19.root'
f_pp = ROOT.TFile(filename_pp, 'read')

filename_AA = 'AA.R04.July18.root'
f_AA = ROOT.TFile(filename_AA, 'read')

for pt_bin in range(n_pt_bins):

    # pp 
    #hname_stat = f'h0_pt{pt_bin}'    # Note this is a TH1
    gname_sys = f'h0_pt{pt_bin}_sys_err'

    g_sys = f_pp.Get(gname_sys)
    x_bins_pp, y_pp, yerr_pp = tgraph_to_numpy(g_sys)


    print()
    print(f'pt bin -- {pt_bin}')
    print(f'  x_bins: {x_bins_pp.tolist()}')
    print(f'  y_pp: {y_pp.tolist()}')
    print(f'  yerr_pp: {yerr_pp.tolist()}')

    # AA
    for cent_bin in range(n_cent_bins):

        # Only consider pt-dependence for 0-10%
        if cent_bin == 0 or pt_bin == 0:

            #hname_stat = f'h{cent_bin}_pt{pt_bin}'    # Note this is a TH1
            gname_sys = f'h{cent_bin}_pt{pt_bin}_sys_err'

            g_sys = f_AA.Get(gname_sys)
            x_bins_AA, y_AA, yerr_AA = tgraph_to_numpy(g_sys)

            if not np.array_equal(x_bins_pp, x_bins_AA):
                print(f'ERROR: x_bins_pp and x_bins_AA are not equal')

            y_ratio = np.divide(y_AA, y_pp)
            yerr_relative_ratio = np.sqrt( np.square( np.divide(yerr_pp, y_pp) ) + np.square( np.divide(yerr_AA, y_AA) ) )
            yerr_ratio = np.multiply(y_ratio, yerr_relative_ratio)

            #print(f'  y_AA -- cent {cent_bin}: {y_AA.tolist()}')
            print(f'  y_ratio -- cent {cent_bin}: {y_ratio.tolist()}')
            print(f'  yerr_ratio -- cent {cent_bin}: {yerr_ratio.tolist()}')