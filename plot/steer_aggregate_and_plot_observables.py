"""
  Macro to steer calculation of observables from histograms produced from a set of runs on XSEDE and uploaded to OSN.

  (See steer_plot_observables.py instead for details on the full workflow of computing observables from final_state_hadrons).

  To run this script:

    - Edit the configuration options at the start of main()

    - Set up environment:
        If downloading from OSN ("download_runinfo" or "download_histograms" True):
            cd STAT-XSEDE-2021/scripts
            python -m venv venv            # only needed the first time
            source venv/bin/activate
            pip install .                  # only needed the first time
        If merging/plotting histograms ("merge_histograms", "aggregate_histograms", or "plot_histograms" True), 
        need ROOT compiled with python, e.g.:
            export ROOTSYS=/path/to/root/root-current
            source $ROOTSYS/bin/thisroot.sh
            (e.g. cd JETSCAPE-analysis && pipenv shell && source init_local.sh)

    - Run script: python plot/steer_aggregate_and_plot_observables.py

------------------------------------------------------------------------

  The workflow is as follows, with each step toggle-able below:

   (1) Download run info for the set of runs on OSN.
       
       This involves two pieces:
         (i) Download runs.yaml for each facility, from STAT-XSEDE-2021/docs/DataManagement
         (ii) Using these runs.yaml, download run_info.yaml for all runs and populate a dictionary with all relevant info for each run
    
   (2) Download histograms for each run.
     
   (3) Merge histograms for each run together.

   (4) Aggregate runs, using the dictionary from step (1).
       
       For each sqrts and design point, merge all histograms into a single one.
       (summed over facilities, run numbers, centralities)

   (5) Plot final observables, and write table for input to Bayesian analysis.

       In the AA case, we plot the AA/pp ratios
       In the pp case, we plot the pp distributions

       Make plots for each design point, as well as a few global QA plots

       We support the same specifications for retrieval of experimental data as described for histogram binning.
       See the STAT_{sqrts}.yaml files for documentation on how to specify the centralities to be looped over.

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import subprocess
import yaml
import pickle
from collections import defaultdict
import pathlib

# ---------------------------------------------------------------
def main():

    #-----------------------------------------------------------------
    # Set which options you want to execute
    download_runinfo = False
    download_histograms = False
    merge_histograms = False
    aggregate_histograms = False
    plot_and_save = True

    # Edit these parameters
    analysis_name = 'Analysis1'
    facilities = ['bridges2', 'expanse']
    stat_xsede_2021_dir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/STAT-XSEDE-2021'
    jetscape_analysis_dir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/JETSCAPE-analysis'
    local_base_outputdir = '/Users/jamesmulligan/JETSCAPE/jetscape-docker/xsede_Analysis1'
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # (1) Download info for all runs from OSN, and create a dictionary with all info needed to download and aggregate histograms
    #-----------------------------------------------------------------
    if download_runinfo or download_histograms or merge_histograms or aggregate_histograms:

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

        histogram_download_location = os.path.join(local_base_outputdir, 'histograms_per_run')
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

                outputdir = os.path.join(local_base_outputdir, f'histograms_per_run/{facility}/{run}')
                inputdir = os.path.join(outputdir, 'histograms')

                if os.path.exists(outputdir):

                    ROOT_filenames = os.listdir(os.path.join(local_base_outputdir, f'histograms_per_run/{facility}/{run}/histograms'))
                    file_list = os.path.join(outputdir, 'files_to_merge.txt')
                    with open(file_list, 'w') as f:
                        for filename in ROOT_filenames:
                            f.write(f'{os.path.join(inputdir, filename)}\n')

                    system = run_dictionary[facility][run]['system']
                    sqrts = run_dictionary[facility][run]['sqrt_s']
                    fname = f'histograms_{system}_{run}_{sqrts}.root'

                    cmd = f'hadd -j 8 -f {os.path.join(outputdir, fname)} @{file_list}'
                    subprocess.run(cmd, check=True, shell=True)
                    os.remove(file_list)

    #-----------------------------------------------------------------
    # Aggregate histograms for runs with a common: (sqrts, system, parametrization_type, design_point_index)
    # We sum over: facilities, run numbers, centralities
    #
    # Write dictionary to design_point_info.pkl that stores relevant info for each design point
    #
    # Note that we store xsec and weight_sum separately for each centrality, 
    # such that we can merge before running plotting script     
    #-----------------------------------------------------------------
    if aggregate_histograms:
        print('Aggregate histograms for runs with common (sqrts, system, parametrization, design_point_index)')
        print()

        # Create a dict that stores list of local paths for each aggregated histogram:
        #   design_point_dictionary[(sqrts, system, parametrization_type, design_point_index)] = ['path/to/run1', 'path/to/run2', ...]
        design_point_dictionary = {}
        for facility in facilities:
            for run in runs[facility]:

                filepath_base = os.path.join(local_base_outputdir, f'histograms_per_run/{facility}/{run}')
                if not list(pathlib.Path(filepath_base).rglob('*.root')):
                    continue

                sqrts = run_dictionary[facility][run]['sqrt_s']
                system = run_dictionary[facility][run]['system']
                if run_dictionary[facility][run]['calculation_type'] == 'jet_energy_loss':
                    parametrization = run_dictionary[facility][run]['parametrization']
                    parametrization_type = parametrization['type']
                    design_point_index = parametrization['design_point_index']
                else:
                    parametrization = None
                    parametrization_type = None
                    design_point_index = None

                design_point_tuple = (sqrts, system, parametrization_type, design_point_index)
                if design_point_tuple not in design_point_dictionary.keys():
                    design_point_dictionary[design_point_tuple] = defaultdict(list)

                filepath = os.path.join(filepath_base, f'histograms_{system}_{run}_{sqrts}.root')
                design_point_dictionary[design_point_tuple]['files'].append(filepath)
                design_point_dictionary[design_point_tuple]['parametrization'] = parametrization

        # Merge each list of histograms together, and write into a new directory structure
        outputdir_base = os.path.join(local_base_outputdir, 'histograms_aggregated')
        for design_point_tuple in design_point_dictionary.keys():
            print(design_point_tuple)
            print()
            sqrts, system, parametrization_type, design_point_index = design_point_tuple

            if parametrization_type:
                fname = f'histograms_design_point_{design_point_index}.root'
                outputdir = os.path.join(outputdir_base, f'{sqrts}_{system}_{parametrization_type}')
            else:
                fname = f'histograms.root'
                outputdir = os.path.join(outputdir_base, f'{sqrts}_{system}')
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)

            cmd = f'hadd -f {os.path.join(outputdir, fname)}'
            for filepath in design_point_dictionary[design_point_tuple]['files']:
                cmd += f' {filepath}'
            subprocess.run(cmd, check=True, shell=True)
            design_point_dictionary[design_point_tuple]['histogram_aggregated'] = os.path.join(outputdir, fname)
            print()

        # Write design_point_dictionary to file
        outfile = os.path.join(outputdir_base, 'design_point_info.pkl')
        with open(outfile, 'wb') as f:
            pickle.dump(design_point_dictionary, f)
        
    #-----------------------------------------------------------------
    # Plot histograms and save output table
    #-----------------------------------------------------------------
    if plot_and_save:

        # Get dictionary containing info for each design point: (sqrts, system, parametrization_type, design_point_index)
        outputdir_base = os.path.join(local_base_outputdir, 'histograms_aggregated')
        outfile = os.path.join(outputdir_base, 'design_point_info.pkl')
        with open(outfile, 'rb') as f:
            design_point_dictionary = pickle.load(f)

        # First plot pp
        for design_point_tuple in design_point_dictionary.keys():
            sqrts, system, parametrization_type, design_point_index = design_point_tuple
            if system == 'pp':
                outputdir = os.path.join(local_base_outputdir, f'plot/{sqrts}_{system}')       
                inputfile = os.path.join(outputdir_base, f'{sqrts}_{system}/histograms.root')       
                cmd = f'python {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                cmd += f' -i {inputfile}'
                cmd += f' -o {outputdir}'
                subprocess.run(cmd, check=True, shell=True)

        # Then plot AA, using appropriate pp reference
        for design_point_tuple in design_point_dictionary.keys():
            sqrts, system, parametrization_type, design_point_index = design_point_tuple
            if system in ['AuAu', 'PbPb']:
                outputdir = os.path.join(local_base_outputdir, f'plot/{sqrts}_{system}_{parametrization_type}')       
                inputfile = os.path.join(outputdir_base, f'{sqrts}_{system}_{parametrization_type}/histograms_design_point_{design_point_index}.root')  
                pp_reference_filename = os.path.join(local_base_outputdir, f'plot/{sqrts}_pp/final_results.root') 
                if os.path.exists(pp_reference_filename):
                    cmd = f'python {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                    cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                    cmd += f' -i {inputfile}'
                    cmd += f' -r {pp_reference_filename}'
                    cmd += f' -o {outputdir}'
                    subprocess.run(cmd, check=True, shell=True)

        # Write PbPb/pp values for each design point to table to be used for Bayesian analysis

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    main()
