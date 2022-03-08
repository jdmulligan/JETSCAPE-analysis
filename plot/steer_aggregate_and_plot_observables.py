"""
  Macro to steer calculation of observables from histograms produced from a set of runs on XSEDE and uploaded to OSN.

  (See steer_plot_observables.py instead for details on the full workflow of computing observables from final_state_hadrons).

  To run this script:

    - Set up environment:
        If downloading from OSN ("download_runinfo" or "download_histograms" True):
            cd STAT-XSEDE-2021/scripts
            python -m venv venv            # only needed the first time
            source venv/bin/activate
            pip install .                  # only needed the first time
        If merging/plotting histograms ("merge_histograms", "aggregate_histograms", "plot_histograms" True), 
        need ROOT compiled with python, e.g.:
            export ROOTSYS=/path/to/root/root-current
            source $ROOTSYS/bin/thisroot.sh
            (e.g. cd JETSCAPE-analysis && pipenv shell && source init_local.sh)
    - Edit the configuration options at the start of main()
    - Run script: python plot/steer_aggregate_and_plot_observables.py

------------------------------------------------------------------------

  The workflow is as follows:

   (1) Download run info for the set of runs on OSN.
       
       This involves two steps:
         (i) Download runs.yaml for each facility, from STAT-XSEDE-2021/docs/DataManagement
         (ii) Using these runs.yaml, download run_info.yaml for all runs and populate a dictionary with all relevant info for each run
    
   (2) Download histograms for each run.
     
   (3) Merge histograms for each run together.

   (4) Aggregate runs, using the dictionary from step (1).
       
       For each sqrts and design point, merge all histograms into a single one.
       (summed over facilities, run numbers, centralities)

   (5) Plot final observables and write table for input to Bayesian analysis.

       In the AA case, we plot the AA/pp ratios
       In the pp case, we plot the pp distributions

       We support the same specifications for retrieval of experimental as described for histogram binning.
       See the STAT_{sqrts}.yaml files for documentation on how to specify the centralities to be looped over.

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import subprocess
import yaml
from collections import defaultdict

# ---------------------------------------------------------------
def main():

    #-----------------------------------------------------------------
    # Set which options you want to execute
    download_runinfo = False
    download_histograms = False
    merge_histograms = True
    aggregate_histograms = False
    plot_histograms = False

    # Edit these parameters
    analysis_name = 'Analysis1'
    facilities = ['bridges2', 'expanse']
    stat_xsede_2021_dir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/STAT-XSEDE-2021'
    local_base_outputdir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/xsede_Analysis1'
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # (1) Download info for all runs from OSN, and create a dictionary with all info needed to download and aggregate histograms
    #-----------------------------------------------------------------
    runs = {}
    run_dictionary = {}

    # (i) Load the runs.yaml for each facility from STAT-XSEDE-2021
    for facility in facilities.copy():
        runs_filename = os.path.join(stat_xsede_2021_dir, f'docs/DataManagement/{analysis_name}/{facility}/runs.yaml')
        if os.path.exists(runs_filename):
            with open(runs_filename, 'r') as f:
                runs[facility] = [name for name in list(yaml.safe_load(f).keys())]     
        else:
            print(f'Warning: runs.yaml not found for {facility}. Removing it from the facilities list.')
            facilities.remove(facility)

    # (ii) Using these runs.yaml, download run_info.yaml for all runs and populate dictionary with relevant info for aggregation
    for facility in facilities:

        run_dictionary[facility] = defaultdict(dict)
        for run in runs[facility]:

            run_info_download_location = os.path.join(local_base_outputdir, 'run_info')
            if download_runinfo:
                download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                cmd = f'python {download_script} -s {facility}/{run}/ -d {run_info_download_location} -f {run}_info.yaml'
                subprocess.run(cmd, check=True, shell=True)

            # Add the run_info block to the run_dictionary
            run_info_file = os.path.join(run_info_download_location, f'{facility}/{run}/{run}_info.yaml')
            if os.path.exists(run_info_file):
                with open(run_info_file, 'r') as f:
                    run_info = yaml.safe_load(f)
                    run_dictionary[facility][run]['calculation_type'] = run_info['calculation_type']
                    run_dictionary[facility][run]['sqrt_s'] = run_info['sqrt_s']
                    run_dictionary[facility][run]['centrality'] = run_info['centrality']
                    if run_info['calculation_type'] == 'jet_energy_loss':
                        if run_info['sqrt_s'] in [200]:
                            run_dictionary[facility][run]['system'] = 'AuAu'
                        elif run_info['sqrt_s'] in [2760, 5020]:
                            run_dictionary[facility][run]['system'] = 'PbPb'
                        run_dictionary[facility][run]['design_point_index'] = run_info['parametrization']['design_point_index']
                        run_dictionary[facility][run]['parametrization'] = run_info['parametrization']
                    else:
                        run_dictionary[facility][run]['system'] = 'pp'
            else:
                if download_runinfo:
                    print(f'Warning: {run}_info.yaml not found on OSN')
            if download_runinfo:
                print()

    # Print what we found
    print('We found the following runs:')
    print()
    for facility in facilities:
        print(f'  {facility}:')
        print()
        for run in runs[facility]:
                print(f'    {run}: {dict(run_dictionary[facility][run])}')
        print()

    #-----------------------------------------------------------------
    # Download histograms for each run
    #-----------------------------------------------------------------
    if download_histograms:
        print('Downloading all histograms...')
        print()

        histogram_download_location = os.path.join(local_base_outputdir, 'histograms')
        for facility in facilities:
            for run in runs[facility]:
                download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                cmd = f'python {download_script} -s {facility}/{run}/ -d {histogram_download_location} -f histograms'
                subprocess.run(cmd, check=True, shell=True)
                print()

        print('Done!')

    #-----------------------------------------------------------------
    # Merge for each run into a single histogram per run
    #-----------------------------------------------------------------
    if merge_histograms:

        for facility in facilities:
            for run in runs[facility]:

                outputdir = os.path.join(local_base_outputdir, f'histograms/{facility}/{run}')
                inputdir = os.path.join(outputdir, 'histograms')

                if os.path.exists(outputdir):

                    ROOT_files = os.listdir(os.path.join(local_base_outputdir, f'histograms/{facility}/{run}/histograms'))
                    system = run_dictionary[facility][run]['system']
                    sqrts = run_dictionary[facility][run]['sqrt_s']
                    fname = f'histograms_{system}_Run{run}_{sqrts}_merged.root'

                    cmd = f'hadd -j 8 -f {os.path.join(outputdir, fname)}'
                    for file in ROOT_files:
                        if '.root' in file:
                            cmd += f' {os.path.join(inputdir, file)}'
                    subprocess.run(cmd, check=True, shell=True)

    #-----------------------------------------------------------------
    # Aggregate histograms for runs with a common sqrts and design point
    # (summed over facilities, run numbers, centralities)
    #
    # Note that we store xsec and weight_sum separately for each centrality, 
    # such that we can merge before running plotting script     
    #-----------------------------------------------------------------
    if aggregate_histograms:
        print('...')

        for run in runlist:
    
            local_run_dir = f'{local_base_dir}/Run{run}'
            inputdir = os.path.join(local_run_dir, 'histograms')
            outputdir = os.path.join(local_run_dir, 'plot')
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)
            
            system = os.listdir(inputdir)[0].split('_')[1]    
            sqrts = os.listdir(inputdir)[0].split('_')[3]
            #centrality = ?
            #design_point_index = ?
            #parameterization_type = ?

            ROOT_files = os.listdir(inputdir)
            fname = f'histograms_{system}_Run{run}_{sqrts}_merged.root'
            cmd = f'hadd -f {os.path.join(outputdir, fname)}'
            for file in ROOT_files:
                if '.root' in file:
                    cmd += f' {os.path.join(inputdir, file)}'
            subprocess.run(cmd, check=True, shell=True)
        
    #-----------------------------------------------------------------
    # Plot histograms
    #-----------------------------------------------------------------
    if plot_histograms:

        # First plot pp
        pp_reference_filenames = {}
        for run in pp_runlist:

            local_run_dir = f'{local_base_dir}/Run{run}'
            inputdir = os.path.join(local_run_dir, 'plot')       
            system = 'pp'   
            sqrts = 5020

            fname = f'histograms_{system}_Run{run}_{sqrts}_merged.root'
            cmd = f'python plot/plot_results_STAT.py -c config/STAT_{sqrts}.yaml -i {inputdir}/{fname}'
            print(cmd)
            subprocess.run(cmd, check=True, shell=True)

            pp_reference_filenames[sqrts] = f'{inputdir}/final_results.root'

        # Then plot AA, using appropriate pp reference
        for run in AA_runlist:

            local_run_dir = f'{local_base_dir}/Run{run}'
            inputdir = os.path.join(local_run_dir, 'plot')       
            system = 'PbPb'   
            sqrts = 5020

            inputdir = os.path.join(local_run_dir, 'plot')
            pp_reference_filename = pp_reference_filenames[sqrts]
            fname = f'histograms_{system}_Run{run}_{sqrts}_merged.root'
            cmd = f'python plot/plot_results_STAT.py -c config/STAT_{sqrts}.yaml -i {inputdir}/{fname} -r {pp_reference_filename}'
            print(cmd)
            subprocess.run(cmd, check=True, shell=True)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    main()
