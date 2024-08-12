"""
  Macro to steer calculation of observables from histograms produced from a set of runs on XSEDE and uploaded to OSN.

  (See steer_plot_observables.py instead for details on the full workflow of computing observables from final_state_hadrons).

  To run this script:

    - Edit the configuration options at the start of main()

    - Set up environment:
        If downloading from OSN ("download_runinfo" or "download_histograms" True):
            cd STAT-XSEDE-2021/scripts
            python3 -m venv venv            # only needed the first time
            source venv/bin/activate
            pip install .                  # only needed the first time
        If merging/plotting histograms ("merge_histograms", "aggregate_histograms", or "plot_histograms" True),
        need ROOT compiled with python, e.g.:
            [enter virtual environment with needed python3 packages, and same python3 version as ROOT build]
            export ROOTSYS=/home/software/users/james/heppy/external/root/root-current
            source $ROOTSYS/bin/thisroot.sh

    - Run script: python3 plot/steer_aggregate_and_plot_observables.py

------------------------------------------------------------------------

  The workflow is as follows, with each step toggle-able below:

   (1) Download run info for the set of runs on OSN.

       This involves two pieces:
         (i) Download runs.yaml for each facility, from STAT-XSEDE-2021/docs/DataManagement
         (ii) Using these runs.yaml, download run_info.yaml for all runs and populate a dictionary with all relevant info for each run

       By default this will only download files that you have not downloaded locally. You can force re-download all with force_download=True.

   (2) Download histograms for each run.

       By default this will only download files that you have not downloaded locally. You can force re-download all with force_download=True.

   (3) Merge histograms for each run together.

   (4) Aggregate runs, using the dictionary from step (1).

       For each (sqrts, system, parametrization_type, design_point_index),
       merge all histograms into a single one, summed over facilities, run numbers, centralities.

   (5) Plot final observables for each design point, and write table for input to Bayesian analysis.

       In the AA case, we plot the AA/pp ratios
       In the pp case, we plot the pp distributions

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import pathlib
import pickle
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

try:
    import ROOT
except ImportError:
    pass

# Suppress performance warning which seems to come from h5py
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# ---------------------------------------------------------------
def main():

    #-----------------------------------------------------------------
    # Set which options you want to execute
    # To edit the set of observables that are plotted, see the self.analysis variable in plot_results_STAT.py
    # Note: all observables that are plotted are saved to hdf5, and are then written to tables
    # Note: plotting script may have warnings about missing histograms, this is usually due to some missing centrality bins for certain design points
    download_runinfo = True
    download_histograms = False
    list_paths_for_selected_design_points = False
    merge_histograms = False
    aggregate_histograms = False
    plot_and_save_histograms = False
    write_tables = False
    plot_global_QA = False

    # TODO keep processing the same run, if we have already processed it, this is useful for debugging
    # TODO update stat-xsede to have updated parsl version?
    # re-analysis parameters
    download_final_state_hadrons = True
    analysis_final_state_hadrons = True # an
    force_reanalysis = True # force to analyse and download again, even if histogram files are present
    delete_final_state_hadrons_after_analysis = True
    # debug options
    ranomize_run_order = True  # if true, runs will be shuffled into random order. this can be useful for benchmarking
    do_debug_same_run = True # if true, only process the same run over and over, useful for debugging



    # Edit these parameters
    stat_xsede_2021_dir = '/software/flo/myJETSCAPE/STAT-XSEDE-2021'
    jetscape_analysis_dir = '/software/flo/myJETSCAPE/JETSCAPE-analysis'
    local_base_outputdir = '/alf/data/flo/jetscape/aggregation_data'
    analyis_container_path = '/alf/data/flo/jetscape/containers/stat_local_gcc_v3.6.sif'
    # analyis_container_path = '/tmp/containers/stat_local_gcc_v3.6.sif'
    force_download = False
    download_threads = 5
    n_cores = 20
    # TODO check skipping logic
    skip_run_binomial = True

    # You may need to edit these for a future analysis -- but can leave as is for now
    analysis_name = 'Analysis3'
    # facilities = ['bridges2', 'expanse'] 
    # facilities = ['test_587cluster']
    facilities = ['bridges2', 'expanse']

    local_analysis_facility = "test_587cluster" # when enabling download_final_state_hadron and subsequent local analysis, specify details about the facility where processing will happen (which is different to the facitliy where these files were downloaded from)


    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # (1) Download info for all runs from OSN, and create a dictionary with all info needed to download and aggregate histograms
    #-----------------------------------------------------------------
    if download_runinfo or download_histograms or list_paths_for_selected_design_points or merge_histograms or aggregate_histograms or download_final_state_hadrons:


        runs = {}
        run_dictionary = {}
        missing_runinfo = defaultdict(list)
        skipped_runs = defaultdict(list) # TODO check why not filled

        # (i) Load the runs.yaml for each facility from STAT-XSEDE-2021
        for facility in facilities.copy():
            runs_filename = os.path.join(stat_xsede_2021_dir, f'docs/DataManagement/{analysis_name}/{facility}/runs.yaml')
            print('Searching for runs.yaml in:', runs_filename)
            if os.path.exists(runs_filename):
                with open(runs_filename, 'r') as f:
                    runs[facility] = [name for name in list(yaml.safe_load(f).keys())]
            else:
                print(f'Warning: runs.yaml not found for {facility}. Removing it from the facilities list.')
                facilities.remove(facility)

        # (ii) Using these runs.yaml, download run_info.yaml for all runs and populate dictionary with relevant info for aggregation
        from js_stat_xsede_steer import file_management
        for facility in facilities:

            # if requested, shuffle runs into random order
            if ranomize_run_order:
                import random
                random.shuffle(runs[facility])

            run_dictionary[facility] = defaultdict(dict)
            max_parallel_downloads = 5
            process_list = []

            # create a multiprocessing pool with 30 workers
            # from multiprocessing import Pool
            # download_pool = Pool(processes=max_parallel_downloads)

            # make a list of tuples that contain source and target
            runinfo_pairs = []

            for i,run in enumerate(runs[facility].copy()):
                
                # if ( i > 2): # todo remove me
                #     break
                run_info_download_location = os.path.join(local_base_outputdir, 'run_info')
                run_info_file = os.path.join(run_info_download_location, f'{facility}/{run}/{run}_info.yaml')
                # print out run_info_file path
                print(run_info_file)

                if download_runinfo:
                    if os.path.exists(run_info_file) and not force_download:
                        print(f'File already exists, will not re-download: {run_info_file} ')
                    else:
                        if force_download:
                            print('Force download enabled -- re-download all runinfo files...')
                        
                        download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                        cmd = f'python3 {download_script} -s {facility}/{run}/ -d {run_info_download_location} -f {run}_info.yaml'

                        # source = f'{facility}/{run}/{run}_info.yaml'
                        source = pathlib.Path(facility) / run / f'{run}_info.yaml'
                        destination = pathlib.Path(run_info_download_location) / f'{facility}/{run}/{run}_info.yaml'
                        # destination = os.path.join(run_info_download_location, f'{facility}/{run}/{run}_info.yaml')
                        destination.parent.mkdir(parents=True, exist_ok=True)


                        # append a tuple with .source and .target to the list
                        # runinfo_pairs.append as FilePair
                        runinfo_pairs.append(file_management.FilePair(source, destination))

                        # add subprocess to pool so that we can run multiple downloads in parallel as specified for number of workers
                        # res = download_pool.apply_async(subprocess.run, args=(cmd,), kwds={'check':True, 'shell':True})

            # check that all processes in pool are done before continuing
            # download_pool.close()
            # download_pool.join()
            failed = file_management.download_from_OSN_pairs(runinfo_pairs, max_parallel_downloads)
            
            if failed:
                print(f'Warning: failed to download run_info.yaml for the following runs: {failed}')
            

            # need to loop over the runs again after processing is finished    
            print('Finished downloading run_info.yaml for all runs, building dictionaries ...')     
            for i,run in enumerate(runs[facility].copy()):

                # should exist now after previous download
                run_info_file = os.path.join(run_info_download_location, f'{facility}/{run}/{run}_info.yaml')
                # Add the run_info block to the run_dictionary
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

                            # if requested, skip runs with binomial parametrization
                            if skip_run_binomial:
                                if run_info['parametrization']['type'] == 'binomial':
                                    # print(f'Skipping run {run} with binomial parametrization...')
                                    skipped_runs[facility].append(run)
                                    runs[facility].remove(run)
                        else:
                            run_dictionary[facility][run]['system'] = 'pp'
                else:
                    print(f'Warning: {run}_info.yaml not found!')
                    runs[facility].remove(run)
                    missing_runinfo[facility].append(run)
                
                # if download_runinfo:
                    # print()

        # Print what we found
        for facility in facilities:
            print(f'{facility}:')
            print()
            print(f'  We found the following runs (N ={len(runs[facility])}):')
            print()
            print(f'    {list(dict(run_dictionary[facility]).keys())}')
            print()
            print(f'Warning: We did NOT find run_info for the following runs:')
            print(f'    {missing_runinfo[facility]}')
            print()
            print()
            print(f'Warning: We skipped the following runs with binomial parametrization:')
            print(f'    {skipped_runs[facility]}')
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

                # If directory exists locally, check number of histograms already downloaded and number on OSN
                histogram_dir = os.path.join(histogram_download_location, f'{facility}/{run}/histograms')
                if os.path.exists(histogram_dir):

                    n_histograms_local = len([h for h in os.listdir(histogram_dir) if '.root' in h])

                    script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/count_files_on_OSN.py')
                    cmd = f'python3 {script} -s {facility}/{run}/ -f histograms'
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                    output = proc.stdout.read()
                    n_histograms_expected = int(output.split()[-1])

                # Download all histograms we are missing (or force download, if requested)
                if os.path.exists(histogram_dir) and n_histograms_expected == n_histograms_local and not force_download:
                    print(f'Histogram dir ({histogram_dir}) already exists and n_histograms_expected ({n_histograms_expected}) == n_histograms_local ({n_histograms_local}), will not re-download')
                    print()
                else:
                    if force_download:
                        print('Force download enabled -- re-download all histograms...')
                    if os.path.exists(histogram_dir) and n_histograms_expected != n_histograms_local:
                        print(f'Histogram dir already exists, but n_histograms_expected ({n_histograms_expected}) does not equal n_histograms_local ({n_histograms_local}), so we will redownload.')

                    download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                    cmd = f'python3 {download_script} -s {facility}/{run}/ -d {histogram_download_location} -f histograms'
                    subprocess.run(cmd, check=True, shell=True)
                    print()

        print('Done!')
    
    if download_final_state_hadrons:
        print('Downloading all final state hadrons...')
        print()

        final_state_hadrons_download_location = os.path.join(local_base_outputdir, 'final_state_hadrons_per_run')
        reananalysis_output_location_base = os.path.join(local_base_outputdir, 'histograms_per_run')
        
        # things for visualization
        total_runs = sum([len(runs[facility]) for facility in facilities])
        run_counter = 0
        for facility in facilities:
            for run in runs[facility]:
                # use for debugging only. Continue to take the same run for testiing
                if do_debug_same_run:
                    run = runs[facility][0]
                run_counter += 1
                run_number = int(run[3:])

                print(f'Processing run {run_number} ({run_counter}/{total_runs}) ...')
                # If directory exists locally, check number of histograms already downloaded and number on OSN
                final_state_hadron_dir = os.path.join(final_state_hadrons_download_location, f'{facility}/{run}')
                reananalysis_output_location = os.path.join(reananalysis_output_location_base, f'{facility}/{run}')
                
                # before we download any analysis output, check if we already have analysis output that we want. If it is already there then skip this run
                if analysis_final_state_hadrons and not force_reanalysis and not do_debug_same_run:
                    histpath = os.path.join(reananalysis_output_location, 'histograms')
                    # before starting analysis check if output directory already exists
                    if os.path.exists(histpath):
                        # check number of output files
                        n_histoutput = len([h for h in os.listdir(histpath) if '.root' in h])
                        n_parquet_local = len([h for h in os.listdir(final_state_hadron_dir) if '.parquet' in h])

                        if n_histoutput == n_parquet_local:
                            print(f'Analysis already done for run {run_number}, I found {n_histoutput} histograms, will not re-analyze ...') 
                            continue

                if os.path.exists(final_state_hadron_dir):

                    n_parquet_local = len([h for h in os.listdir(final_state_hadron_dir) if '.parquet' in h])

                    script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/count_files_on_OSN.py')
                    cmd = f'python3 {script} -s {facility}/{run}/ -f final_state_hadrons'
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                    output = proc.stdout.read()
                    n_parquet_expected = int(output.split()[-1])

                # Download all files we are missing (or force download, if requested)
                if os.path.exists(final_state_hadron_dir) and n_parquet_expected == n_parquet_local and not force_download:
                    print(f'Hadron file dir ({final_state_hadron_dir}) already exists and n_parquet_expected ({n_parquet_expected}) == n_parquet_local ({n_parquet_local}), will not re-download')
                    print()
                else:
                    if force_download or do_debug_same_run:
                        print('Force download enabled -- re-download all parquet files...')
                    if os.path.exists(final_state_hadron_dir) and n_parquet_expected != n_parquet_local:
                        print(f'Histogram dir already exists, but n_parquet_expected ({n_parquet_expected}) does not equal n_parquet_local ({n_parquet_local}), so we will redownload.')

                    download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                    cmd = f'python3 {download_script} -s {facility}/{run}/ -d {final_state_hadrons_download_location} -f final_state_hadrons -c {download_threads}'
                    subprocess.run(cmd, check=True, shell=True)
                    print()
                
                # run analyis script for that run
                # this will go here, we will also need to grab an analysis config, maybe possible rom xstat repo directly? maybe one can also just create is herej
                # we are still in run loop, from here on onwards pseudo code since i do not know how to make the codes talk to each other
                # from submit.py run function run_analysis_standalone
                # def run_analysis_standalone(
                # run_number: int, -> implemented!
                # run_info_dir: Path, -> implemented!
                # final_state_hadrons_dir: Path, -> implemented
                # output_dir: Path, -> implemented
                # facility_name: str, -> implemented
                # container_path: Path, -> implemented
                # facility_config: parsl_configs.Facility, ->TODO
                # local_jetscape_analysis_override_path: Optional[Path] = None, -> implemented


                # instead of config, just pass facility name (which should be enough to singularity binary and options)
                # run analhysis script, do not know yet how to load
                # make sure to install it

                # run is string RunXXXXX, extract only runnumber by removing Run from string
                if analysis_final_state_hadrons:
                    from js_stat_xsede_steer import submit
                    ntasks = 80
                    max_n_slurm_jobs = 80

                    submit.run_analysis_standalone(run_number,run_info_download_location, final_state_hadrons_download_location, reananalysis_output_location, local_analysis_facility,facility, analyis_container_path,ntasks, max_n_slurm_jobs, jetscape_analysis_dir)

                # TODO, take care of deleting the final_state_hadrons files locally for this run after analysis is done. Maybe do this as an option, unsure what exactly I will do
                # get parquet files for final state hadrons final_state_hadrons_dir / facility / RunXXXXX
                # run_directory = final_state_hadrons_dir / facility_name / f"Run{run_number:04}"

                if delete_final_state_hadrons_after_analysis:
                    # safely delete final_state_hadron_dir, make sure the path contains final_state_hadrons_per_run so we don't by accident delete "/" or "/home" .... 

                    if os.path.exists(final_state_hadron_dir) and "final_state_hadrons_per_run" in final_state_hadron_dir:
                        # remove all files in the directory that contain final_state_hadrons and parquet in name
                        # note that we don't ask for endging parquet since sometimes incomplete files have weird ending
                        print(f'Deleting all parquet files in {final_state_hadron_dir} as requested...')
                        for f in os.listdir(final_state_hadron_dir):
                            if "final_state_hadrons" in f and "parquet" in f:
                                os.remove(os.path.join(final_state_hadron_dir, f))
                        # now we can delete the directory
                        os.rmdir(final_state_hadron_dir)

                # TODO upload analysis results / histos to OSN in the future
                


        print('Done!')
    

    #-----------------------------------------------------------------
    # Merge for each run into a single histogram per run
    #-----------------------------------------------------------------
    if merge_histograms:

        for facility in facilities:
            run_counter = 0
            for run in runs[facility]:
                run_counter += 1

                # TODO make counter "Processing run X out of Y at facility Z"
                print(f'Processing run {run_counter} out of {len(runs[facility])} at facility {facility}...')

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

                    cmd = f'hadd -j {n_cores} -f {os.path.join(outputdir, fname)} @{file_list}'
                    subprocess.run(cmd, check=True, shell=True)
                    os.remove(file_list)

    #-----------------------------------------------------------------
    # List histogram paths for selected design points.
    # Used to extract unmerged histograms for separate studies.
    #-----------------------------------------------------------------
    if list_paths_for_selected_design_points:
        # These options are only put here since they're quite niche.
        selected_parametrization = "exponential"
        selected_design_point_indices = list(range(40))
        selected_centrality = [0, 10]
        # Using a relative directory here is useful since we will want to tar these files up
        # while keeping a directory structure
        relative_dir = Path(local_base_outputdir)

        # Setup
        # NOTE: We'll store a list of paths per design index because we have multiple sqrt_s per design point
        design_point_index_to_path = defaultdict(list)

        for facility in facilities:
            for run in runs[facility]:

                filepath_base = os.path.join(local_base_outputdir, f'histograms_per_run/{facility}/{run}')
                if not list(pathlib.Path(filepath_base).rglob('*.root')):
                    continue

                sqrts = run_dictionary[facility][run]['sqrt_s']
                system = run_dictionary[facility][run]['system']
                if run_dictionary[facility][run]['calculation_type'] == 'jet_energy_loss':
                    # AA case
                    centrality = run_dictionary[facility][run]['centrality']
                    parametrization = run_dictionary[facility][run]['parametrization']
                    parametrization_type = parametrization['type']
                    design_point_index = parametrization['design_point_index']

                    # Apply selection
                    if not (
                        centrality == selected_centrality
                        and design_point_index in selected_design_point_indices
                        and parametrization_type == selected_parametrization
                    ):
                        continue
                else:
                    # pp case
                    # We define -1 as the convention
                    design_point_index = -1
                    continue

                design_point_index_to_path[design_point_index].append(Path(filepath_base))

        # Sort the output
        # NOTE: We don't sort by sqrt_s because we don't have that info. It's not so important in any case.
        design_point_index_to_path = dict(sorted(design_point_index_to_path.items()))
        output_files = []
        for paths in design_point_index_to_path.values():
            for path in paths:
                # We want two outputs:
                # 1. The merged histogram, if available.
                # NOTE: We need to search for the histogram before we convert it to the relative path!
                output_files.extend([p.relative_to(relative_dir) for p in path.glob("histograms_*.root")])
                # 2. The unmerged histograms directory.
                output_files.append(path.relative_to(relative_dir) / "histograms")

        # And print it for the user to utilize as desired
        print("We found the following paths for the selected design points:\n")
        print()
        print(design_point_index_to_path)
        print()
        print("  List of files:")
        print("    " + " ".join([str(s) for s in output_files]))


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
    # Plot histograms and save to ROOT files
    #-----------------------------------------------------------------
    if plot_and_save_histograms:

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
                cmd = f'python3 {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                cmd += f' -i {inputfile}'
                cmd += f' -o {outputdir}'
                subprocess.run(cmd, check=True, shell=True)

        # Then plot AA, using appropriate pp reference, and construct AA/pp ratios
        process_list = []
        for design_point_tuple in design_point_dictionary.keys():
            sqrts, system, parametrization_type, design_point_index = design_point_tuple
            if system in ['AuAu', 'PbPb']:
                outputdir = os.path.join(local_base_outputdir, f'plot/{sqrts}_{system}_{parametrization_type}/{design_point_index}')
                inputfile = os.path.join(outputdir_base, f'{sqrts}_{system}_{parametrization_type}/histograms_design_point_{design_point_index}.root')
                pp_reference_filename = os.path.join(local_base_outputdir, f'plot/{sqrts}_pp/final_results.root')
                if os.path.exists(pp_reference_filename):
                    cmd = f'python3 {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                    cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                    cmd += f' -i {inputfile}'
                    cmd += f' -r {pp_reference_filename}'
                    cmd += f' -o {outputdir}'

                    # Execute in parallel
                    # Simple & quick implementation: once max_processes have been launched, wait for them to finish before continuing
                    process = subprocess.Popen(cmd, shell=True)
                    process_list.append(process)
                    if len(process_list) > n_cores-1:
                        for subproc in process_list:
                            subproc.wait()
                        process_list = []

    #-----------------------------------------------------------------
    # Convert histograms to data tables and save
    #-----------------------------------------------------------------
    if write_tables:

        # Construct table of model predictions for all design points, as input to Bayesian analysis
        # We can easily adapt this to the format v1.0 specified here, although we may want to update a bit: https://www.evernote.com/l/ACWFCWrEcPxHPJ3_P0zUT74nuasCoL_DBmY
        plot_dir = os.path.join(local_base_outputdir, 'plot')
        table_base_dir = os.path.join(local_base_outputdir, 'tables')
        prediction_table_dir = os.path.join(table_base_dir, 'Prediction')
        data_table_dir = os.path.join(table_base_dir, 'Data')
        design_table_dir = os.path.join(table_base_dir, 'Design')

        if not os.path.exists(prediction_table_dir):
            os.makedirs(prediction_table_dir)
        if not os.path.exists(data_table_dir):
            os.makedirs(data_table_dir)
        if not os.path.exists(design_table_dir):
            os.makedirs(design_table_dir)

        # We will write out design point files as well
        design_point_dir = os.path.join(local_base_outputdir, 'histograms_aggregated')
        outfile = os.path.join(design_point_dir, 'design_point_info.pkl')
        with open(outfile, 'rb') as f:
            design_point_dictionary = pickle.load(f)
        design_df = {}

        # Loop through each directory corresponding to a given (sqrts, parameterization)
        for label in os.listdir(plot_dir):
            if 'AuAu' in label or 'PbPb' in label:
                sqrts, system, parameterization  = label.split('_')

                output_dict = defaultdict()
                output_dict['values'] = {}
                output_dict['errors'] = {}
                output_dict['bin_edges'] = {}
                output_dict['observable_label'] = {}

                # Loop through design points and observables, and construct a dataframe of predictions for each observable
                #   columns=[design1, design2, ...]
                #   rows=[bin1, bin2, ...]
                label_dir = os.path.join(plot_dir, label)
                print(f'Constructing prediction dataframes for {label}...')
                for design_point_index in os.listdir(label_dir):
                    #if design_point_index != '0':
                    #    continue
                    if not os.path.isfile(os.path.join(label_dir, design_point_index)) and design_point_index != 'Data':
                        #print(f'  design_point_index: {design_point_index}')

                        final_result_h5 = os.path.join(f'{label_dir}/{design_point_index}', 'final_results.h5')
                        if not os.path.exists(final_result_h5):
                            continue

                        with h5py.File(final_result_h5, 'r') as hf:
                            for key in list(hf.keys()):

                                # Use a separate dataframe for values, errors, bin_edges
                                if 'values' in key:
                                    type = 'values'
                                elif 'errors' in key:
                                    type = 'errors'
                                elif 'bin_edges' in key:
                                    type = 'bin_edges'
                                else:
                                    sys.exit(f'Unexpected key: {key}')

                                # Get observable label for bookkeeping
                                observable_labels = [s.replace('h_','',1).replace('.pdf','') for s in os.listdir(os.path.join(label_dir, design_point_index)) if 'pdf' in s]
                                output_dict['observable_label'][key] = None
                                for s in observable_labels:
                                    if s in key:
                                        output_dict['observable_label'][key] = s

                                # Put design point info into dataframe, with design point index as index of dataframe
                                if key not in output_dict[type]:
                                    output_dict[type][key] = pd.DataFrame()

                                output_dict[type][key][f'design_point{design_point_index}'] = hf[key][:]

                                design_point_key = (int(sqrts), system, parameterization, int(design_point_index))
                                parameterization_values = pd.DataFrame(data=[design_point_dictionary[design_point_key]['parametrization']['parametrization_values']],
                                                                       index=[int(design_point_index)])
                                if parameterization not in design_df.keys():
                                    design_df[parameterization] = parameterization_values
                                else:
                                    if int(design_point_index) not in design_df[parameterization].index:
                                        design_df[parameterization] = pd.concat([design_df[parameterization], parameterization_values])
                print(f'Done constructing prediction dataframes for {label}.')
                print()

                # Write Prediction and Data dataframes to txt
                print(f'Writing prediction tables for {label}...')
                for type in output_dict.keys():
                    if type in ['observable_label', 'bin_edges']:
                        continue

                    for key,df in output_dict[type].items():
                        key_items = key.split('_')
                        #if 'values' in key:
                        #    print()
                        #    print(key_items)

                        # Parse observable-specific names -- there are a few different cases depending on the observable class
                        # We uniformize the structure as: f'{parameterization}__{sqrts}__{system}__{observable_category}__{observable}__{subobservable}__{centrality[0]}-{centrality[1]}'
                        # For example: binomial__5020__PbPb__inclusive_chjet__pt_alice__R0.4__0-10
                        # This will allow us to easily parse the observables in a uniform way, and also access info in the STAT observable config files
                        # Note that the experiment can always be accessed as observable.split('_')[-1]
                        # The subobservable can include: jet radius, grooming condition, pt bin index

                        # Hadron observables -- after hole subtraction
                        if 'hadron' in key and 'unsubtracted' not in key:
                            observable_category = key_items[2]
                            observable = f'{key_items[3]}_{key_items[4]}_{key_items[5]}'
                            subobservable = ''
                            centrality = [''.join(filter(str.isdigit, s)) for s in key_items[6].split(',')]

                        # Jet observables -- with negative_recombiner
                        # There are several different subobservable patterns that we need to parse
                        elif 'negative_recombiner' in key:
                            if 'zcut' in key:
                                observable_category = f'{key_items[4]}_{key_items[5]}'
                                observable = f'{key_items[6]}_{key_items[7]}'
                                subobservable = f'{key_items[8]}_{key_items[9]}_{key_items[10]}'
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[11].split(',')]
                            elif 'charge' in key:
                                observable_category = f'{key_items[4]}_{key_items[5]}'
                                observable = f'{key_items[6]}_{key_items[7]}'
                                subobservable = f'{key_items[8]}_{key_items[9]}'
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[10].split(',')]
                            elif 'dijet' in key:
                                observable_category = key_items[4]
                                observable = f'{key_items[5]}_{key_items[6]}'
                                subobservable = f'{key_items[7]}_{key_items[9]}'
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[8].split(',')]
                            elif '_y_' in key:
                                observable_category = f'{key_items[4]}_{key_items[5]}'
                                observable = f'{key_items[6]}_{key_items[7]}_{key_items[8]}'
                                subobservable = f'{key_items[9]}_{key_items[11]}'
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[10].split(',')]
                            elif key_items[-2] in [f'pt{i}' for i in range(10)]:
                                observable_category = f'{key_items[4]}_{key_items[5]}'
                                observable = f'{key_items[6]}_{key_items[7]}'
                                subobservable = f'{key_items[8]}_{key_items[10]}'
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[9].split(',')]
                            else:
                                observable_category = f'{key_items[4]}_{key_items[5]}'
                                observable = f'{key_items[6]}_{key_items[7]}'
                                subobservable = key_items[8]
                                centrality = [''.join(filter(str.isdigit, s)) for s in key_items[9].split(',')]

                        # Skip other observables (hadrons without subtraction, jets with other subtraction schemes)
                        else:
                            if 'shower_recoil' not in key and 'constituent_subtraction' not in key and 'unsubtracted' not in key:
                                print(f'Unexpected key: {key}')
                            continue

                        observable_name_data = f'{sqrts}__{system}__{observable_category}__{observable}__{subobservable}__{centrality[0]}-{centrality[1]}'
                        observable_name_prediction = f'{parameterization}__{observable_name_data}'

                        if 'values' in key:
                            print(f'  {observable_name_prediction}')

                        # Sort columns
                        df_prediction = output_dict[type][key]
                        df_prediction = df_prediction.reindex(sorted(df_prediction.columns, key=lambda x: float(x[12:])), axis=1)
                        columns = list(df_prediction.columns)

                        # Get experimental data
                        observable_label = output_dict['observable_label'][key]
                        data_dir = os.path.join(plot_dir, f'{label}/Data')
                        filename = f'Data_{observable_label}.dat'
                        data_file = os.path.join(data_dir, filename)
                        if os.path.exists(data_file):
                            data = np.loadtxt(data_file, ndmin=2)

                        if df_prediction.to_numpy().shape[0] != data.shape[0]:
                            print(f'Mismatch of number of bins: prediction ({df_prediction.to_numpy().shape[0]}) vs. data ({data.shape[0]})')

                        # Remove rows with leading zeros (corresponding to bins below the min_pt cut)
                        n_zero_rows = 0
                        for row in df_prediction.to_numpy():
                            if np.all(row == 0):
                                n_zero_rows += 1
                            else:
                                break
                        df_prediction = df_prediction.iloc[n_zero_rows:, :]
                        data = data[n_zero_rows:, :]

                        if df_prediction.to_numpy().shape[0] != data.shape[0]:
                            print(f'Mismatch of number of bins after removing zeros: prediction ({df_prediction.to_numpy().shape[0]}) vs. data ({data.shape[0]})')
                            print()

                        # Write Prediction.dat and Data.dat files
                        filename = os.path.join(prediction_table_dir, f'Prediction__{observable_name_prediction}__{type}.dat')
                        design_point_file = f'Design_{parameterization}.dat'
                        header = f'Version 2.0\nData Data_{observable_name_prediction}.dat\nDesign {design_point_file}\n'
                        header += ' '.join(columns)
                        np.savetxt(filename, df_prediction.values, header=header)

                        filename = os.path.join(data_table_dir, f'Data__{observable_name_data}.dat')
                        header = f'Version 1.1\n'
                        header += 'Label xmin xmax y y_err'
                        np.savetxt(filename, data, header=header)

                print(f'Done writing prediction tables for {label}.')
                print()

        # Write out the Design.dat files
        print('Writing design point tables...')
        for parameterization in design_df.keys():

            # Sort according to index
            df = design_df[parameterization].sort_index()

            # Rename columns
            df.rename(columns={'t_start': 'Tau0', 'alpha_s': 'AlphaS', 'q_switch': 'Q0', 'A': 'C1',  'B': 'C2'}, inplace=True)
            if parameterization == 'binomial':
                df.rename(columns={'C': 'A', 'D': 'B'}, inplace=True)
            elif parameterization == 'exponential':
                df.rename(columns={'C': 'C3'}, inplace=True)

            # Reorder columns to match previous convention
            if parameterization == 'binomial':
                ordered_columns = ['AlphaS', 'Q0', 'C1', 'C2', 'Tau0', 'A', 'B']
            elif parameterization == 'exponential':
                ordered_columns = ['AlphaS', 'Q0', 'C1', 'C2', 'Tau0', 'C3']
            df = df[ordered_columns]

            # Write
            filename = os.path.join(design_table_dir, f'Design__{parameterization}.dat')
            header = f'Version 1.0\n'
            header += f'- Design points for {parameterization} PDF\n'
            parameters = ' '.join(df.keys())
            header += f'Parameter {parameters}\n'
            header += f'- Parameter AlphaS: Linear [0.1, 0.5]\n'
            header += f'- Parameter Q0: Linear [1, 10]\n'
            header += f'- Parameter C1: Log [0.006737946999085467, 10]\n'
            header += f'- Parameter C2: Log [0.006737946999085467, 10]\n'
            header += f'- Parameter Tau0: Linear [0.0, 1.5]\n'
            if parameterization == 'binomial':
                header += f'- Parameter A: Linear [-10, 100]\n'
                header += f'- Parameter B: Linear [-10, 100]\n'
            elif parameterization == 'exponential':
                header += f'- Parameter C3: Log [0.049787068367863944, 100]\n'
            indices = ' '.join([str(i) for i in df.index])
            header += f'Design point indices (row index): {indices}'
            np.savetxt(filename, df.values, header=header)

        print('Done!')

    #-----------------------------------------------------------------
    # Plot some QA over all aggregated runs
    #-----------------------------------------------------------------
    if plot_global_QA:

        histograms_aggregated_dir = os.path.join(local_base_outputdir, 'histograms_aggregated')
        plot_dir = os.path.join(local_base_outputdir, 'plot')
        global_qa_dir = os.path.join(local_base_outputdir, 'global_qa')
        if not os.path.exists(global_qa_dir):
            os.makedirs(global_qa_dir)

        #----------------------
        # Plot n_events
        # Make a 2d plot from a 2d numpy array: n_generated/n_target for (sqrts_paramterization_centrality, design_point_index)
        n_design_points = 230
        n_design_point_max = {'exponential': 230, 'binomial': 180}
        n_systems = 12
        shape = (n_systems, n_design_points)
        n_events = np.zeros(shape)
        n_target = {'200': 500000, '2760': 442000, '5020': 1700000}

        # Keep track of which design points need to be re-run
        rerun_threshold = 0.75
        rerun_dict = defaultdict(list)
        missing_dict = defaultdict(list)

        # Loop through each directory corresponding to a given (sqrts, parameterization)
        system_index = 0
        system_labels = []
        for sqrts_parameterization_label in os.listdir(plot_dir):
            if 'AuAu' in sqrts_parameterization_label or 'PbPb' in sqrts_parameterization_label:
                sqrts, system, parameterization  = sqrts_parameterization_label.split('_')
                qa_plot_dir = os.path.join(plot_dir, sqrts_parameterization_label)

                for design_point_index in range(n_design_points):

                    if design_point_index >= n_design_point_max[parameterization]:
                        continue

                    if not os.path.exists(os.path.join(qa_plot_dir, str(design_point_index))):
                        missing_dict[f'{sqrts_parameterization_label}_0-10'].append(int(design_point_index))
                        missing_dict[f'{sqrts_parameterization_label}_10-50'].append(int(design_point_index))
                        continue

                    fname = f'{sqrts_parameterization_label}/histograms_design_point_{design_point_index}.root'
                    f = ROOT.TFile(os.path.join(histograms_aggregated_dir, fname), 'read')
                    h = f.Get('h_centrality_generated')
                    ratio_0_10 = h.Integral(h.GetXaxis().FindBin(0+0.5), h.GetXaxis().FindBin(10-0.5)) / n_target[sqrts]
                    ratio_10_50 = h.Integral(h.GetXaxis().FindBin(10+0.5), h.GetXaxis().FindBin(50-0.5)) / n_target[sqrts]

                    n_events[system_index, int(design_point_index)] = ratio_0_10
                    n_events[system_index+1, int(design_point_index)] = ratio_10_50

                    if 0.01 < ratio_0_10 < rerun_threshold:
                        rerun_dict[f'{sqrts_parameterization_label}_0-10'].append(int(design_point_index))
                    if 0.01 < ratio_10_50 < rerun_threshold:
                        rerun_dict[f'{sqrts_parameterization_label}_10-50'].append(int(design_point_index))
                    if np.isclose(ratio_0_10, 0.):
                        missing_dict[f'{sqrts_parameterization_label}_0-10'].append(int(design_point_index))
                    if np.isclose(ratio_10_50, 0.):
                        missing_dict[f'{sqrts_parameterization_label}_10-50'].append(int(design_point_index))

                system_labels.append(f'{sqrts_parameterization_label}_0-10')
                system_labels.append(f'{sqrts_parameterization_label}_10-50')
                system_index += 2

        # Order the systems
        ordered_indices = [3, 2, 11, 10, 9, 8, 5, 4, 7, 6, 1, 0]
        n_events[:] = n_events[ordered_indices,:]
        system_labels = [system_labels[i] for i in ordered_indices]

        # Plot
        fig, ax = plt.subplots()
        fig.suptitle(r'Number of events: $N_{gen} / N_{target}$', fontsize=16)
        c = ax.imshow(n_events, cmap='jet', aspect='auto', vmin=0., vmax=2., interpolation='nearest')
        fig.colorbar(c)
        ax.set_xlabel('Design point index', size=12)
        system_ticks = np.linspace(0, n_systems-1, n_systems)
        plt.yticks(system_ticks, system_labels, size=6)
        outfilename = os.path.join(global_qa_dir, f'n_events.pdf')
        plt.tight_layout()
        plt.savefig(outfilename)
        print('Plotted n_events.')

        # Print which design points needs to be rerun
        print()
        print(f'The following design points had fewer than {rerun_threshold}x the target statistics')
        for key,val in rerun_dict.items():
            print(f'  {key}: {sorted(val)}')
        print()
        print(f'The following design points have not yet been uploaded')
        for key,val in missing_dict.items():
            print(f'  {key}: {sorted(val)}')
        print()

        #----------------------
        # Plot statistical uncertainty for each bin of each observable (for a given design point)
        table_dir = os.path.join(local_base_outputdir, 'tables')
        prediction_dir = os.path.join(table_dir, 'Prediction')
        data_dir = os.path.join(table_dir, 'Data')

        parameterizations = ['exponential', 'binomial']
        n_bins = 31

        # Construct 3d arrray: relative statistical uncertainty for (design_point_index, observable, bin)
        for parameterization in parameterizations:

            n_observables = int(len([x for x in os.listdir(prediction_dir) if 'Prediction' in x and 'values' in x and parameterization in x]))
            shape = (n_bins, n_design_point_max[parameterization], n_observables)
            relative_uncertainty_prediction = np.zeros(shape)
            relative_uncertainty_ratio_to_data = np.zeros(shape)

            observable_labels = []
            i_observable = 0
            for file_prediction in os.listdir(prediction_dir):
                #file_prediction_keys = file_prediction.split('__')
                if 'Prediction' in file_prediction and 'values' in file_prediction and parameterization in file_prediction:
                    observable = file_prediction[12:-12].replace(f'{parameterization}__', '')
                    file_prediction_errors = file_prediction.replace('values', 'errors')

                    # Get predictions and compute relative uncertainty -- zero pad to fixed size
                    prediction_values = np.loadtxt(os.path.join(prediction_dir, file_prediction), ndmin=2)
                    prediction_errors = np.loadtxt(os.path.join(prediction_dir, file_prediction_errors), ndmin=2)
                    if 0 in prediction_values:
                        print(f'WARNING: {file_prediction} has value=0 at design points {np.where(prediction_values == 0)[1]}')

                    relative_uncertainty_prediction_unpadded = np.divide(prediction_errors, prediction_values)
                    if relative_uncertainty_prediction_unpadded.shape[0] > n_bins:
                        sys.exit(f'Set n_bins to be {relative_uncertainty_prediction_unpadded.shape[0]} or larger (due to {observable})')
                    observable_shape = relative_uncertainty_prediction_unpadded.shape
                    n_pads_x = shape[0] - observable_shape[0]
                    n_pads_y = shape[1] - observable_shape[1]
                    relative_uncertainty_prediction_padded = np.pad(relative_uncertainty_prediction_unpadded, ((0,n_pads_x), (0,n_pads_y)))
                    relative_uncertainty_prediction[:,:,i_observable] = relative_uncertainty_prediction_padded

                    # Get data and compute relative uncertainty -- zero pad to fixed size
                    file_data = f'Data__{observable}.dat'
                    data = np.loadtxt(os.path.join(data_dir, file_data), ndmin=2)
                    data_values = data[:,2]
                    data_errors = data[:,3]
                    if 0 in data_values:
                        print(f'WARNING: {file_data} has value=0 at design points {np.where(prediction_values == 0)[1]}')

                    relative_uncertainty_data_unpadded = np.divide(data_errors, data_values)

                    # Compute ratio of prediction uncertainty to data uncertainty -- zero pad to fixed size
                    if prediction_values.shape[0] != data_values.shape[0]:
                        sys.exit(f'({observable_shape}) has different shape than Data ({relative_uncertainty_data_unpadded.shape})')

                    relative_uncertainty_ratio_to_data_unpadded = np.divide(relative_uncertainty_prediction_unpadded,
                                                                            relative_uncertainty_data_unpadded[:,None])
                    relative_uncertainty_ratio_to_data_padded = np.pad(relative_uncertainty_ratio_to_data_unpadded, ((0,n_pads_x), (0,n_pads_y)))
                    relative_uncertainty_ratio_to_data[:,:,i_observable] = relative_uncertainty_ratio_to_data_padded

                    observable_labels.append(observable)
                    i_observable += 1

            # TODO: Do something about bins that have value=0 and cause division problem, leading to empty entry

            # Order the observables by sqrts, and then alphabetically (i.e. by observable and centrality)
            ordered_labels = np.sort(observable_labels).tolist()
            ordered_indices = np.argsort(observable_labels).tolist()
            ordered_indices_200 = [ordered_indices[i] for i in range(n_observables) if '200' in ordered_labels[i]]
            ordered_indices_2700 = [ordered_indices[i] for i in range(n_observables) if '2760' in ordered_labels[i]]
            ordered_indices_5020 = [ordered_indices[i] for i in range(n_observables) if '5020' in ordered_labels[i]]
            ordered_indices = ordered_indices_5020 + ordered_indices_2700 + ordered_indices_200

            relative_uncertainty_prediction[:] = relative_uncertainty_prediction[:,:,ordered_indices]
            relative_uncertainty_ratio_to_data[:] = relative_uncertainty_ratio_to_data[:,:,ordered_indices]
            observable_labels_ordered = [observable_labels[i] for i in ordered_indices]
            print()
            print(f'(n_bins, n_design_points, n_observables_total) = {relative_uncertainty_prediction.shape}')

            # Group observables and make a plot for each group
            groups = ['hadron__pt_', 'jet__pt_', 'Dpt', 'Dz', '']
            group_indices_total = []
            for group in groups:
                if group: # Plot all observables containing a matching string
                    group_indices = [i for i,label in enumerate(observable_labels_ordered) if group in label]
                    group_indices_total += group_indices
                else: # Plot everything that remains
                    group_indices = [i for i,_ in enumerate(observable_labels_ordered) if i not in group_indices_total]
                    group = 'other'
                n_observables_group = len(group_indices)
                print(f'n_observables ({group}) = {n_observables_group}')

                group_mask = np.array([i in group_indices for i,_ in enumerate(observable_labels_ordered)])
                observable_labels_ordered_group = [observable_labels_ordered[i] for i in group_indices]

                relative_uncertainty_prediction_group = relative_uncertainty_prediction[:,:,group_mask]
                relative_uncertainty_ratio_to_data_group = relative_uncertainty_ratio_to_data[:,:,group_mask]

                # Plot relative uncertainty of prediction
                fig = plt.figure(figsize=[12, 15])
                ax = plt.axes()
                fig.suptitle(f'% statistical uncertainty on prediction -- mean', fontsize=24)
                matrix = np.transpose(np.mean(relative_uncertainty_prediction_group, axis=1))
                matrix_masked = np.ma.masked_where((matrix < 1e-8), matrix)
                c = ax.imshow(matrix_masked, cmap='jet', aspect='auto', vmin=0., vmax=0.2, interpolation='nearest')
                fig.colorbar(c)
                ax.set_xlabel('Observable bin', size=16)
                bin_ticks = range(n_bins)
                plt.xticks(bin_ticks, bin_ticks, size=10)
                observable_ticks = np.linspace(0, n_observables_group-1, n_observables_group)
                plt.yticks(observable_ticks, observable_labels_ordered_group, size=10)
                outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_prediction_{parameterization}_{group}_mean.pdf')
                plt.tight_layout()
                plt.savefig(outfilename)
                plt.close()

                # Plot ratio of relative uncertainty of prediction to that in data
                fig = plt.figure(figsize=[12, 15])
                ax = plt.axes()
                fig.suptitle(f'(% statistical uncertainty on prediction) / (% total uncertainty on data) -- mean', fontsize=18)
                matrix = np.transpose(np.mean(relative_uncertainty_ratio_to_data_group, axis=1))
                matrix_masked = np.ma.masked_where((matrix < 1e-8), matrix)
                c = ax.imshow(matrix_masked, cmap='jet', aspect='auto', vmin=0., vmax=5., interpolation='nearest')
                fig.colorbar(c)
                ax.set_xlabel('Observable bin', size=16)
                bin_ticks = range(n_bins)
                plt.xticks(bin_ticks, bin_ticks, size=10)
                observable_ticks = np.linspace(0, n_observables_group-1, n_observables_group)
                plt.yticks(observable_ticks, observable_labels_ordered_group, size=10)
                outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_ratio_{parameterization}_{group}_mean_0_5.pdf')
                plt.tight_layout()
                plt.savefig(outfilename)
                plt.close()

                # Repeat with different z axis range
                fig = plt.figure(figsize=[12, 15])
                ax = plt.axes()
                fig.suptitle(f'(% statistical uncertainty on prediction) / (% total uncertainty on data) -- mean', fontsize=18)
                matrix = np.transpose(np.mean(relative_uncertainty_ratio_to_data_group, axis=1))
                matrix_masked = np.ma.masked_where((matrix < 1e-8), matrix)
                c = ax.imshow(matrix_masked, cmap='jet', aspect='auto', vmin=0., vmax=1., interpolation='nearest')
                fig.colorbar(c)
                ax.set_xlabel('Observable bin', size=16)
                bin_ticks = range(n_bins)
                plt.xticks(bin_ticks, bin_ticks, size=10)
                observable_ticks = np.linspace(0, n_observables_group-1, n_observables_group)
                plt.yticks(observable_ticks, observable_labels_ordered_group, size=10)
                outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_ratio_{parameterization}_{group}_mean_0_1.pdf')
                plt.tight_layout()
                plt.savefig(outfilename)
                plt.close()

                # Plot for each design point
                plot_each_design_point = False
                for design_point_index in range(n_design_points):

                    if plot_each_design_point or design_point_index in [0]:

                        # Plot relative uncertainty of prediction
                        fig = plt.figure(figsize=[12, 15])
                        ax = plt.axes()
                        fig.suptitle(f'% statistical uncertainty on prediction -- design point {design_point_index}', fontsize=24)
                        matrix = np.transpose(relative_uncertainty_prediction[:,design_point_index,group_mask])
                        matrix_masked = np.ma.masked_where((matrix < 1e-8), matrix)
                        c = ax.imshow(matrix_masked, cmap='jet', aspect='auto', vmin=0., vmax=0.2, interpolation='nearest')
                        fig.colorbar(c)
                        ax.set_xlabel('Observable bin', size=16)
                        bin_ticks = range(n_bins)
                        plt.xticks(bin_ticks, bin_ticks, size=10)
                        observable_ticks = np.linspace(0, n_observables_group-1, n_observables_group)
                        plt.yticks(observable_ticks, observable_labels_ordered_group, size=10)
                        outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_prediction_{parameterization}_{group}_design_point{design_point_index}.pdf')
                        plt.tight_layout()
                        plt.savefig(outfilename)
                        plt.close()

                        # Plot ratio of relative uncertainty of prediction to that in data
                        fig = plt.figure(figsize=[12, 15])
                        ax = plt.axes()
                        fig.suptitle(f'(% statistical uncertainty on prediction) / (% total uncertainty on data) -- design point {design_point_index}', fontsize=18)
                        matrix = np.transpose(relative_uncertainty_ratio_to_data[:,design_point_index,group_mask])
                        matrix_masked = np.ma.masked_where((matrix < 1e-8), matrix)
                        c = ax.imshow(matrix_masked, cmap='jet', aspect='auto', vmin=0., vmax=1.0, interpolation='nearest')
                        fig.colorbar(c)
                        ax.set_xlabel('Observable bin', size=16)
                        bin_ticks = range(n_bins)
                        plt.xticks(bin_ticks, bin_ticks, size=10)
                        observable_ticks = np.linspace(0, n_observables_group-1, n_observables_group)
                        plt.yticks(observable_ticks, observable_labels_ordered_group, size=10)
                        outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_ratio_{parameterization}_{group}_design_point{design_point_index}.pdf')
                        plt.tight_layout()
                        plt.savefig(outfilename)
                        plt.close()

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    main()
