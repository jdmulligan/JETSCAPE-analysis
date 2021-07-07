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

    sqrts_list = ['200', '2760', '5020']
    construct_observables = False
    construct_histograms = False
    merge_histograms = False
    plot_histograms = True

    final_state_hadron_dirs = {}
    for sqrts in sqrts_list:
        final_state_hadron_dirs[sqrts] = f'/Users/jamesmulligan/JETSCAPE/jetscape-docker/JETSCAPE-output/STAT/pp_{sqrts}_final_state_hadrons'
    
    for sqrts,dir in final_state_hadron_dirs.items():
    
        #-----------------------------------------------------------------
        # Loop through final_state_hadron files, and construct observables
        if construct_observables:

            # Clear output directory
            outputdir = dir.replace('final_state_hadrons', 'observables')
            if os.path.exists(outputdir):
                shutil.rmtree(outputdir)

            for file in os.listdir(dir):
                cmd = f'python jetscape_analysis/analysis/analyze_events_STAT.py -c config/STAT_{sqrts}.yaml -i {dir}/{file} -o {outputdir}'
                print(cmd)
                subprocess.run(cmd, check=True, shell=True)
          
        #-----------------------------------------------------------------
        # Loop through observable files, and construct histograms
        if construct_histograms:
        
            inputdir = dir.replace('final_state_hadrons', 'observables')
            outputdir = dir.replace('final_state_hadrons', 'histograms')

            for file in os.listdir(inputdir):
                if 'observables' in file:
                    cmd = f'python plot/histogram_results_STAT.py -c config/STAT_{sqrts}.yaml -i {inputdir}/{file} -o {outputdir}'
                    print(cmd)
                    subprocess.run(cmd, check=True, shell=True)
    
        #-----------------------------------------------------------------
        # Merge histograms
        if merge_histograms:
        
            # Merge
            inputdir = dir.replace('final_state_hadrons', 'histograms')
            outputdir = dir.replace('final_state_hadrons', 'plot')
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
        
            inputdir = dir.replace('final_state_hadrons', 'plot')
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
