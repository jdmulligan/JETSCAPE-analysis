# JETSCAPE-analysis

This repository contains basic tools to generate and analyze Jetscape events in HepMC3 format.

It is written entirely in python -- leveraging c++ underneath where necessary -- no compilation necessary!

## (1) Generating events

The script `generate/generate_jetscape_events.py` generates JETSCAPE events, 
including automated machinery to launch a set of pt-hat bins.

### Pre-requisites

To generate JETSCAPE events, you must first build the [JETSCAPE package](https://github.com/JETSCAPE/JETSCAPE) itself.

We recommend to follow the [JETSCAPE Docker Instructions](https://github.com/JETSCAPE/JETSCAPE/tree/master/docker) to do so. 

Assuming you have a Jetscape docker installation according to the above instructions 
(with a shared folder located at `~/jetscape-docker`, containing the Jetscape repository at `~/jetscape-docker/JETSCAPE`), 
you should do (from outside the docker container):

```
cd ~/jetscape-docker/
git clone git@github.com:jdmulligan/JETSCAPE-analysis.git
```

You should then enter the docker container as specified in the above instructions.
The generation script should then be run from inside the JETSCAPE docker container:

```
python generate_jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml -o /my/outputdir
```

where `jetscapeAnalysisConfig.yaml` should be edited to specify the pt-hat bins and JETSCAPE XML configuration paths, 
and `outputdir` specifies where the JETSCAPE output files will be written.
Note that the machinery here only modifies the pt-hat bins in the JETSCAPE XML configuration -- all other settings should
be set manually (which modules to include, output format type, etc.).

That's it! The script will write a separate sub-directory with JETSCAPE events for each pt-hat bin. 

## (2) Analyzing events

The script `analysis/analyze_jetscape_events.py` analyzes JETSCAPE events, producing an output ROOT file.
It also contains machinery to aggregate the results from the set of pt-hat bins, and plot the analysis results.

### Pre-requisites

Once the JETSCAPE events are generated, we no longer rely on the JETSCAPE package nor its docker container. 
Instead, we analyze the events (jet-finding, writing histograms, etc.) using a local python environment. 
For jet-finding, we rely on the package [heppy](https://github.com/matplo/heppy) which wraps fastjet and fastjet-contribs in python

#### One-time setup

We recommend to use [pipenv(https://github.com/pypa/pipenv) to manage your python environment:

```
cd JETSCAPE-analysis
pipenv --three
pipenv shell
pipenv install pyhepmc_ng pyyaml numpy tqdm ROOT
```

Install `heppy` wherever you desire, and load its modules:

```
cd <my-heppy-location>
git clone git@github.com:matplo/heppy.git
cd heppy
./scripts/setup.sh --buildext --root
module use <my-heppy-location>/modules
module load heppy/main_python
```

#### Workflow

Once you have done the one-time setup, your general workflow is the following:

```
cd JETSCAPE-analysis
pipenv shell
module use <my-heppy-location>/modules
module load heppy/main_python
```

And then run the script:

```
python generate_jetscape_events.py -c ../config/jetscapeAnalysisConfig.yaml -o /my/outputdir
```

where `/my/outputdir` is the directory containing the generated JETSCAPE events.

