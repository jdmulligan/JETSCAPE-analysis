#!/usr/bin/env python3

# General
import os

# Analysis
import numpy as np
import matplotlib.pyplot as plt

def plot_pthat():

    #power_list = [3,4,5,6,7,8]
    power_list = [8]
    pt_ref = 10
    for pow in power_list: 

        output_dir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/JETSCAPE-output/STAT_v2/pt_weight_test'
        filename_weighted = os.path.join(output_dir, f'test_out_final_state_hadrons_pow{pow}_ref{pt_ref}_pthat5-50_prec15-2.dat')
        filename_unweighted = os.path.join(output_dir, f'test_out_final_state_hadrons_powOFF_100K.dat')

        # Get pt_hat array from each file
        pthat_alphaN, weights_alphaN = pthat_array(filename_weighted)
        pthat_alpha0, weights_alpha0 = pthat_array(filename_unweighted)

        # Plot pt-hat array
        plot(pthat_alphaN, weights_alphaN, pthat_alpha0, weights_alpha0, pow, pt_ref, output_dir)

def plot(pthat_alphaN, weights_alphaN, pthat_alpha0, weights_alpha0, pow, pt_ref, output_dir):

    sum_weights_alphaN = np.sum(weights_alphaN)
    sum_weights_alpha0 = np.sum(weights_alpha0)
    print(f'sum_weights_alpha{pow}: {sum_weights_alphaN}')
    print(f'sum_weights_alpha0: {sum_weights_alpha0}')

    fig, (ax1, ax2) = plt.subplots(nrows=2)
    plt.xlabel(r'$\hat{p}_{T}$', fontsize=14)
    alpha = 0.7
    linewidth = 2
    bins = np.linspace(0, 50, 25)
    if pt_ref == 10:
        ax1.set_ylim([1e-5, np.power(10,4+(pow-3))])
    else:
        ax1.set_ylim([1e-10, np.power(10,4+(pow-3)-10.)])
    #ax2.set_ylim([5e-1, 5e4])
    ax2.set_ylim([5*np.power(10, -1*(pt_ref/10.)), 5*np.power(10, 5-(pt_ref/10.))])

    # Plot 1/N_event dN/dx (equal to dsigma/dx, up to sigma_pthat_min factor)
    weights1 = [1./sum_weights_alphaN for _ in weights_alphaN]
    n_weighted,_,_ = ax1.hist(pthat_alphaN,
                            bins,
                            histtype='step',
                            label=rf'$\alpha={pow}$',
                            weights=weights1,
                            linewidth=linewidth,
                            linestyle='-',
                            alpha=alpha)

    weights2 = [w/sum_weights_alphaN for w in weights_alphaN]
    n_weighted_scaled,_,_ = ax1.hist(pthat_alphaN,
                            bins,
                            histtype='step',
                            label=rf'$\alpha={pow}$, weighted',
                            weights=weights2,
                            linewidth=linewidth,
                            linestyle='--',
                            alpha=alpha)

    weights3 = [1./sum_weights_alpha0 for _ in weights_alpha0]
    n_unweighted,_,_ = ax1.hist(pthat_alpha0,
                            bins,
                            histtype='step',
                            label = r'$\alpha=0$',
                            weights=weights3,
                            linewidth=linewidth,
                            linestyle='-',
                            alpha=alpha)
    legend = ax1.legend(loc='best', fontsize=12, frameon=False)
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\frac{1}{\sigma_{\hat{p}_T}} \frac{d\sigma}{d\hat{p}_T}$', fontsize=12)

    # Plot ratio of weighted to unweighted
    ratio = n_weighted/n_unweighted
    ratio[np.isnan(ratio)] = 0
    x = (bins[1:] + bins[:-1]) / 2
    ax2.plot(x, ratio, label=r'{}'.format('$\\frac{{\\alpha={}}} {{\\alpha=0}}$'.format(pow)), linewidth=linewidth, alpha=alpha)

    # Plot inverse of weight
    ax2.plot(x, np.power(x/pt_ref, pow), label=r'{}'.format('$\\left( \\frac{{\\hat{{p}}_T}} {{ p_{{T}}^{{ref}} }} \\right)^{}$'.format(pow)), linewidth=linewidth, alpha=alpha)

    # Plot ratio of weighted_scaled to unweighted
    ratio = n_weighted_scaled/n_unweighted
    ratio[np.isnan(ratio)] = 0
    ax2.plot(x, ratio, label=r'{}'.format('$\\frac{{\\alpha={},weighted}} {{\\alpha=0}}$'.format(pow)), linestyle='--', linewidth=linewidth, alpha=alpha)

    ax2.set_yscale('log')
    ax2.set_ylabel('ratio', fontsize=12)
    legend = ax2.legend(loc='best', fontsize=12, frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join('/Users/jamesmulligan/JETSCAPE/jetscape-docker/JETSCAPE-output/STAT_v2/pt_weight_test', f'pt_hat_pow{pow}_ref{pt_ref}.pdf'))
    plt.close()

    # Also plot histogram of weights
    plt.xscale('log')
    min = np.amin(weights_alphaN)
    max = np.amax(weights_alphaN)
    print(f'min: {min}, max: {max}')
    plt.hist(weights_alphaN,
            bins=np.logspace(np.log10(min),np.log10(max), 100),
            histtype='step',
            density=True,
            label = rf'$\alpha={pow}$',
            linewidth=linewidth,
            linestyle='-',
            alpha=alpha,
            log=True)
    plt.savefig(os.path.join(output_dir, f'weights_pow{pow}_ref{pt_ref}.pdf'))
    plt.close()

def pthat_array(filename):

    pthat_list = []
    weight_list = []

    file = open(filename, 'r')
    for line in file.readlines():
        row = line.split()
        if 'pt_hat' in row:
            pthat_list.append(float(row[-1]))
            weight_list.append(float(row[4]))

    return np.array(pthat_list), np.array(weight_list)

##################################################################
if __name__ == "__main__":
    plot_pthat()