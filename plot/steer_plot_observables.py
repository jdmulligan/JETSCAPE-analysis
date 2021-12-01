"""
  macro to steer calculation and plotting of observables from final_state_hadrons
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import subprocess
import shutil

# ---------------------------------------------------------------
def main():

    sqrts = 5020
    final_state_hadron_dir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/stampede/Run0001'

    construct_observables = False
    construct_histograms = False
    merge_histograms = False
    plot_histograms = True
    
    #-----------------------------------------------------------------
    # Loop through final_state_hadron files, and construct observables
    if construct_observables:

        # Clear output directory
        observables_dir = os.path.join(final_state_hadron_dir, 'observables')
        if os.path.exists(observables_dir):
            shutil.rmtree(observables_dir)

        for file in os.listdir(final_state_hadron_dir):
            if 'final_state_hadrons' in file:
                cmd = f'python jetscape_analysis/analysis/analyze_events_STAT.py -c config/STAT_{sqrts}.yaml -i {final_state_hadron_dir}/{file} -o {observables_dir}'
                print(cmd)
                subprocess.run(cmd, check=True, shell=True)
        
    #-----------------------------------------------------------------
    # Loop through observable files, and construct histograms
    if construct_histograms:
    
        inputdir = os.path.join(final_state_hadron_dir, 'observables')
        outputdir = os.path.join(final_state_hadron_dir, 'histograms')

        for file in os.listdir(inputdir):
            if 'observables' in file:
                cmd = f'python plot/histogram_results_STAT.py -c config/STAT_{sqrts}.yaml -i {inputdir}/{file} -o {outputdir}'
                print(cmd)
                subprocess.run(cmd, check=True, shell=True)

    #-----------------------------------------------------------------
    # Merge histograms
    if merge_histograms:
    
        # Merge
        inputdir = os.path.join(final_state_hadron_dir, 'histograms')
        outputdir = os.path.join(final_state_hadron_dir, 'plot')
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        
        ROOT_files = os.listdir(inputdir)
        fname = f'histograms_{sqrts}_merged.root'
        cmd = f'hadd -f {os.path.join(outputdir, fname)}'
        for file in ROOT_files:
            if '.root' in file:
                cmd += f' {os.path.join(inputdir, file)}'
        subprocess.run(cmd, check=True, shell=True)
        
    #-----------------------------------------------------------------
    # Plot histograms
    if plot_histograms:
    
        inputdir = os.path.join(final_state_hadron_dir, 'plot')
        fname = f'histograms_{sqrts}_merged.root'
        cmd = f'python plot/plot_results_STAT.py -c config/STAT_{sqrts}.yaml -i {inputdir}/{fname}'
        print(cmd)
        subprocess.run(cmd, check=True, shell=True)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    main()
