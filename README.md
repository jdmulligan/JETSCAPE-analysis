# JETSCAPE-analysis

This repository contains basic tools to generate and analyze Jetscape events in HepMC3 and/or Ascii format.

It is written entirely in python – leveraging c++ underneath where necessary – no compilation necessary!

## (1) Generating events

The script `jetscape_analysis/generate/jetscape_events.py` generates JETSCAPE events, including automated machinery to
launch a set of pt-hat bins and optionally scan over any additional parameter(s).

### Pre-requisites

To generate JETSCAPE events, you must first build the [JETSCAPE package](https://github.com/JETSCAPE/JETSCAPE) itself. 
We recommend to follow the [JETSCAPE Docker
Instructions](https://github.com/JETSCAPE/JETSCAPE/tree/master/docker) to do so.

Assuming you have a Jetscape docker installation according to the above instructions
(with a shared folder located at `~/jetscape-docker`, containing the Jetscape repository at `~/jetscape-docker/JETSCAPE`),
you should do (from outside the docker container):

```
cd ~/jetscape-docker/
git clone git@github.com:jdmulligan/JETSCAPE-analysis.git
```

You should then enter the docker container as specified in the above instructions.

### Generate events

The generation script should then be run from inside the JETSCAPE docker container:

```
python jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml -o /home/jetscape-user/JETSCAPE-analysis-output
```

where 
- `-c` specifies a configuration file that should be edited to specify the pt-hat bins and JETSCAPE XML configuration paths,
- `-o` specifies a location where the JETSCAPE output files will be written.

Note that the machinery here only modifies the pt-hat bins and (optionally) other parameter values in the JETSCAPE XML configuration -- but does not allow to change which modules are present.

That's it! The script will write a separate sub-directory with JETSCAPE events for each pt-hat bin.

## (2) Analyzing events

The script `jetscape_analysis/analysis/analyze_events.py` analyzes JETSCAPE events, producing an output ROOT file.
It also contains machinery to aggregate the results from the set of pt-hat bins, and plot the analysis results.

### Pre-requisites

Once the JETSCAPE events are generated, we no longer rely on the JETSCAPE package,
but rather we analyze the events (jet-finding, writing histograms, etc.) using a python environment.
A preconfigured environment is available in the JETSCAPE docker container -- or it can be installed manually.
For jet-finding, we rely on the package [heppy](https://github.com/matplo/heppy) which wraps fastjet 
and fastjet-contribs in python.

#### Docker installation (recommended)

Assuming you have set up Docker according to Step (1), 
enter the container and run an initialization script:

```
cd /home/jetscape-user/JETSCAPE-analysis
source init.sh
```

That's it! Now you can proceed to analyze events.

#### Manual installation

We recommend to use a virtual environment such as [pipenv](https://github.com/pypa/pipenv) to
manage your python environment:

```
cd /home/jetscape-user/JETSCAPE-analysis
pipenv --three
pipenv install pyhepmc_ng pyyaml numpy tqdm
```

Install `heppy` wherever you desire:

```
cd <my-heppy-location>
git clone git@github.com:matplo/heppy.git
cd heppy
./external/build.sh
```

Then load the heppy module:

```
cd /home/jetscape-user/JETSCAPE-analysis
pipenv shell
module use <my-heppy-location>/modules
module load heppy/1.0
```

### Analyze events

Simply run the script:

```
cd /home/jetscape-user/JETSCAPE-analysis/jetscape_analysis/analysis
python analyze_events.py -c ../../config/jetscapeAnalysisConfig.yaml -i /home/jetscape-user/JETSCAPE-analysis-output -o /my/outputdir
```

where 
- `-c` specifies a configuration file that should be edited to specify the pt-hat bins and analysis parameters,
- `-i` specifies  is the directory containing the generated JETSCAPE events,
- `-o` specifies a location where the analysis output will be written.

---------------------------------------------------------------------

### Setup `poetry`

This section is only relevant if you want to package up the code yourself.

We use `poetry` to manage packaging up this code. You need to setup poetry once globally.

You must use version 1.0 or later (currently in 1.0.0b3 as of Nov 2019). Follow the [installation
instructions](https://poetry.eustace.io/docs/#installation). If the python version that you would like to use
doesn't automatically show up in your path (for example, if you are a macOS user and don't use pyenv), you can
set the default python for poetry with (python 3.7 in this example):

```bash
$ poetry env use 3.7
# Check that it worked with:
$ poetry env info
```

You can run commands within the poetry virtual environment using `poetry run <command>`. If you want to load
the virtual environment directly, you can try `poetry shell`, which may work (it doesn't for me), or use
`source "$(dirname $(poetry run which python))/activate"`, as suggested
[here](https://github.com/sdispater/poetry/issues/571#issuecomment-443595960). You can also alias this command
for convenience. As long as it is run within the repository, it will always load the right virtual
environment.

Install this package using `poetry install` from the repository root.

#### Pre-commit checks

To setup checks that run on every commit, run

```bash
$ poetry run pre-commit install
```

Now, each commit will be checked on the users' machine.

