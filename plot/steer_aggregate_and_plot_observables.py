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
            [enter virtual environment with needed python packages, and same python version as ROOT build]
            export ROOTSYS=/home/software/users/james/heppy/external/root/root-current
            source $ROOTSYS/bin/thisroot.sh

    - Run script: python plot/steer_aggregate_and_plot_observables.py

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
import sys
import subprocess
import yaml
import pickle
from collections import defaultdict
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

try:
    import ROOT
except ImportError:
    pass

# ---------------------------------------------------------------
def main():

    #-----------------------------------------------------------------
    # Set which options you want to execute
    download_runinfo = False
    download_histograms = False
    merge_histograms = False
    aggregate_histograms = False
    plot_and_save_histograms = False
    write_tables = False
    plot_global_QA = True

    # Edit these parameters
    stat_xsede_2021_dir = '/home/james/jetscape-docker/STAT-XSEDE-2021'
    jetscape_analysis_dir = '/home/james/jetscape-docker/JETSCAPE-analysis'
    local_base_outputdir = '/rstorage/jetscape/STAT-Bayesian/Analysis1/20220601'
    force_download = False
    n_cores = 20

    # You may need to edit these for a future analysis -- but can leave as is for now
    analysis_name = 'Analysis1'
    facilities = ['bridges2', 'expanse']
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    # (1) Download info for all runs from OSN, and create a dictionary with all info needed to download and aggregate histograms
    #-----------------------------------------------------------------
    if download_runinfo or download_histograms or merge_histograms or aggregate_histograms:

        runs = {}
        run_dictionary = {}
        missing_runinfo = defaultdict(list)

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
            for run in runs[facility].copy():

                run_info_download_location = os.path.join(local_base_outputdir, 'run_info')
                run_info_file = os.path.join(run_info_download_location, f'{facility}/{run}/{run}_info.yaml')

                if download_runinfo:
                    if os.path.exists(run_info_file) and not force_download:
                        print(f'File already exists, will not re-download: {run_info_file} ')
                    else:
                        if force_download:
                            print('Force download enabled -- re-download all runinfo files...')

                        download_script = os.path.join(stat_xsede_2021_dir, f'scripts/js_stat_xsede_steer/download_from_OSN.py')
                        cmd = f'python {download_script} -s {facility}/{run}/ -d {run_info_download_location} -f {run}_info.yaml'
                        subprocess.run(cmd, check=True, shell=True)

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
                        else:
                            run_dictionary[facility][run]['system'] = 'pp'
                else:
                    if download_runinfo:
                        print(f'Warning: {run}_info.yaml not found on OSN')
                    runs[facility].remove(run)
                    missing_runinfo[facility].append(run)

                if download_runinfo:
                    print()

        # Print what we found
        for facility in facilities:
            print(f'{facility}:')
            print()
            print('  We found the following runs:')
            print()
            print(f'    {list(dict(run_dictionary[facility]).keys())}')
            print()
            print(f'Warning: We did NOT find run_info for the following runs:')
            print(f'    {missing_runinfo[facility]}')
            print()
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

                histogram_dir = os.path.join(histogram_download_location, f'{facility}/{run}/histograms')
                if os.path.exists(histogram_dir) and not force_download:
                    print(f'Histogram dir already exists, will not re-download: {histogram_dir} ')
                else:
                    if force_download:
                        print('Force download enabled -- re-download all histograms...')

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

                    cmd = f'hadd -j {n_cores} -f {os.path.join(outputdir, fname)} @{file_list}'
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
                cmd = f'python {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                cmd += f' -i {inputfile}'
                cmd += f' -o {outputdir}'
                subprocess.run(cmd, check=True, shell=True)

        # Then plot AA, using appropriate pp reference, and construct dataframe of PbPb/pp values to be used for Bayesian analysis
        max_processes = 20
        process_list = []
        for design_point_tuple in design_point_dictionary.keys():
            sqrts, system, parametrization_type, design_point_index = design_point_tuple
            if system in ['AuAu', 'PbPb']:
                outputdir = os.path.join(local_base_outputdir, f'plot/{sqrts}_{system}_{parametrization_type}/{design_point_index}')       
                inputfile = os.path.join(outputdir_base, f'{sqrts}_{system}_{parametrization_type}/histograms_design_point_{design_point_index}.root')  
                pp_reference_filename = os.path.join(local_base_outputdir, f'plot/{sqrts}_pp/final_results.root') 
                if os.path.exists(pp_reference_filename):
                    cmd = f'python {jetscape_analysis_dir}/plot/plot_results_STAT.py'
                    cmd += f' -c {jetscape_analysis_dir}/config/STAT_{sqrts}.yaml'
                    cmd += f' -i {inputfile}'
                    cmd += f' -r {pp_reference_filename}'
                    cmd += f' -o {outputdir}'

                    # Execute in parallel
                    # Simple & quick implementation:once max_processes have been launched, wait for them to finish before continuing
                    process = subprocess.Popen(cmd, shell=True)
                    process_list.append(process)
                    if len(process_list) > max_processes:
                        for subproc in process_list:
                            subproc.wait()
                        process_list = []

    #-----------------------------------------------------------------
    # Convert histograms to data tables and save
    #-----------------------------------------------------------------
    if write_tables:

        # Construct table of model predictions for all design points, as input to Bayesian analysis
        # We can easily adapt this to the format v1.0 specified here, although we may want to update a bit: https://www.evernote.com/l/ACWFCWrEcPxHPJ3_P0zUT74nuasCoL_DBmY
        print()
        print('Write predictions to table...')
        table_base_dir = os.path.join(local_base_outputdir, 'tables')
        plot_dir = os.path.join(local_base_outputdir, 'plot')

        if not os.path.exists(table_base_dir):
            os.makedirs(table_base_dir)

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

                # Loop through design points and observables, and construct a dataframe for each observable
                #   columns=[design1, design2, ...]
                #   rows=[bin1, bin2, ...]
                label_dir = os.path.join(plot_dir, label)
                for design_point_index in os.listdir(label_dir):
                    if not os.path.isfile(os.path.join(label_dir, design_point_index)) and design_point_index != 'Data':

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

                # Write Prediction and Data dataframes to txt
                # TODO: For now we use negative recombiner for jet observables
                for type in output_dict.keys():
                    if type == 'observable_label':
                        continue

                    for key,df in output_dict[type].items():

                        key_items = key.split('_')
                        if 'hadron' in key and 'unsubtracted' not in key:
                            experiment = key_items[5]
                            observable = f'{key_items[2]}_{key_items[3]}_{key_items[4]}'
                            centrality = [''.join(filter(str.isdigit, s)) for s in key_items[6].split(',')]
                            observable_name = f'{experiment}_{system}{sqrts}{parameterization}_{observable}_{centrality[0]}to{centrality[1]}'
                        elif 'negative_recombiner' in key and '_y_' not in key:
                            experiment = key_items[7]
                            observable = f'{key_items[4]}_{key_items[5]}_{key_items[6]}_{key_items[8]}'
                            centrality = [''.join(filter(str.isdigit, s)) for s in key_items[9].split(',')]
                            observable_name = f'{experiment}_{system}{sqrts}{parameterization}_{observable}_{centrality[0]}to{centrality[1]}'
                        else:
                            continue

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
                        filename = os.path.join(table_base_dir, f'Prediction_{observable_name}_{type}.dat')
                        design_point_file = f'Design_{parameterization}.dat'
                        header = f'Version 2.0\nData Data_{observable_name}.dat\nDesign {design_point_file}\n'
                        header += ' '.join(columns)
                        np.savetxt(filename, df_prediction.values, header=header)

                        #filename = os.path.join(table_base_dir, f'Data_{observable_name}.dat')
                        #header = f'Version 1.1\n'
                        #header += 'Label xmin xmax y y_err'
                        #np.savetxt(filename, data, header=header)

        # Write out the Design.dat files
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
            filename = os.path.join(table_base_dir, f'Design_{parameterization}.dat')
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
        n_design_points = 140
        n_systems = 12
        shape = (n_systems, n_design_points)
        n_events = np.zeros(shape)
        n_target = {'200': 500000, '2760': 442000, '5020': 1700000}

        # Keep track of which design points need to be re-run
        rerun_threshold = 0.75
        rerun_dict = defaultdict(list)

        # Loop through each directory corresponding to a given (sqrts, parameterization)
        system_index = 0
        system_labels = []
        for sqrts_parameterization_label in os.listdir(plot_dir):
            if 'AuAu' in sqrts_parameterization_label or 'PbPb' in sqrts_parameterization_label:
                sqrts, system, parameterization  = sqrts_parameterization_label.split('_')
                qa_plot_dir = os.path.join(plot_dir, sqrts_parameterization_label)

                for design_point_index in os.listdir(qa_plot_dir):
                    if design_point_index == 'Data':
                        continue
            
                    fname = f'{sqrts_parameterization_label}/histograms_design_point_{design_point_index}.root'
                    f = ROOT.TFile(os.path.join(histograms_aggregated_dir, fname), 'read')
                    h = f.Get('h_centrality_generated')
                    ratio_0_10 = h.Integral(0, 10) / n_target[sqrts]
                    ratio_10_50 = h.Integral(10, 50) / n_target[sqrts]

                    n_events[system_index, int(design_point_index)] = ratio_0_10
                    n_events[system_index+1, int(design_point_index)] = ratio_10_50

                    if ratio_0_10 < rerun_threshold:
                        rerun_dict[f'{sqrts_parameterization_label}_0-10'].append(design_point_index)
                    if ratio_10_50 < rerun_threshold:
                        rerun_dict[f'{sqrts_parameterization_label}_10-50'].append(design_point_index)

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

        #----------------------
        # Plot statistical uncertainty for each bin of each observable (for a given design point)

        table_dir = os.path.join(local_base_outputdir, 'tables')

        parameterizations = ['exponential', 'binomial']
        n_bins = 21
        n_observables = int(len([x for x in os.listdir(table_dir) if 'Prediction' in x and 'values' in x]) / len(parameterizations))
        shape = (n_bins, n_design_points, n_observables)
        relative_uncertainty = np.zeros(shape)

        # Construct 3d arrray: relative statistical uncertainty for (design_point_index, observable, bin)
        for parameterization in parameterizations:
            observable_labels = []
            i_observable = 0
            for file in os.listdir(table_dir):
                if 'Prediction' in file and 'values' in file and parameterization in file:
                    observable = file[11:-11]
                    file_errors = file.replace('values', 'errors') 

                    values = np.loadtxt(os.path.join(table_dir, file), ndmin=2)
                    errors = np.loadtxt(os.path.join(table_dir, file_errors), ndmin=2)
                    relative_uncertainty_unpadded = np.divide(errors, values)

                    # Zero pad to fixed size
                    if relative_uncertainty_unpadded.shape[0] > n_bins:
                        sys.exit(f'Set n_bins to be {relative_uncertainty_unpadded.shape[0]} or larger (due to {observable})')

                    observable_shape = relative_uncertainty_unpadded.shape
                    n_pads_x = shape[0] - observable_shape[0]
                    n_pads_y = shape[1] - observable_shape[1]
                    relative_uncertainty_padded = np.pad(relative_uncertainty_unpadded, ((0,n_pads_x), (0,n_pads_y)))

                    relative_uncertainty[:,:,i_observable] = relative_uncertainty_padded
                    observable_labels.append(observable.replace(parameterization, ''))
                    i_observable += 1

            # Order the observables by sqrts, and then alphabetically (i.e. by observable and centrality)
            ordered_labels = np.sort(observable_labels).tolist()
            ordered_indices = np.argsort(observable_labels).tolist()
            ordered_indices_200 = [ordered_indices[i] for i in range(n_observables) if '200' in ordered_labels[i]]
            ordered_indices_2700 = [ordered_indices[i] for i in range(n_observables) if '2760' in ordered_labels[i]]
            ordered_indices_5020 = [ordered_indices[i] for i in range(n_observables) if '5020' in ordered_labels[i]]
            ordered_indices = ordered_indices_5020 + ordered_indices_2700 + ordered_indices_200

            relative_uncertainty[:] = relative_uncertainty[:,:,ordered_indices]
            observable_labels_ordered = [observable_labels[i] for i in ordered_indices]

            # Plot for each design point
            # TODO: Normalize by uncertainty in experimental measurement
            # TODO: plot average over all design points? (need to exclude points we haven't run yet)
            for design_point_index in range(n_design_points):
                fig = plt.figure(figsize=[12, 15])
                ax = plt.axes()
                fig.suptitle(f'Relative statistical uncertainty, design point {design_point_index}', fontsize=24)
                c = ax.imshow(np.transpose(relative_uncertainty[:,design_point_index,:]), cmap='jet', aspect='auto', vmin=0., vmax=0.2, interpolation='nearest')
                fig.colorbar(c)
                ax.set_xlabel('Observable bin', size=16)
                bin_ticks = range(n_bins)
                plt.xticks(bin_ticks, bin_ticks, size=10)
                observable_ticks = np.linspace(0, n_observables-1, n_observables)
                plt.yticks(observable_ticks, observable_labels_ordered, size=10)
                outfilename = os.path.join(global_qa_dir, f'stat_uncertainty_{parameterization}_design_point{design_point_index}.pdf')
                plt.tight_layout()
                plt.savefig(outfilename)
                plt.close()

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print()
    main()
