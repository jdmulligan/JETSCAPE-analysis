# JETSCAPE-analysis

This repository contains basic tools to generate and analyze Jetscape events in HepMC3 format.

It is written entirely in python -- leveraging c++ underneath where necessary -- no compilation necessary!

## (1) Generating events

The script `generate/jetscape_events.py` generates JETSCAPE events, including automated machinery to
launch a set of pt-hat bins.

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
The generation script should then be run from inside the JETSCAPE docker container:

```
python jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml -o /my/outputdir
```

where `jetscapeAnalysisConfig.yaml` should be edited to specify the pt-hat bins and JETSCAPE XML configuration paths,
and `outputdir` specifies where the JETSCAPE output files will be written.
Note that the machinery here only modifies the pt-hat bins in the JETSCAPE XML configuration -- all other settings should
be set manually (which modules to include, output format type, etc.).

That's it! The script will write a separate sub-directory with JETSCAPE events for each pt-hat bin.

## (2) Analyzing events

The script `analysis/analyze_events.py` analyzes JETSCAPE events, producing an output ROOT file.
It also contains machinery to aggregate the results from the set of pt-hat bins, and plot the analysis results.

### Pre-requisites

Once the JETSCAPE events are generated, we no longer rely on the JETSCAPE package nor its docker container.
Instead, we analyze the events (jet-finding, writing histograms, etc.) using a local python environment.
For jet-finding, we rely on the package [heppy](https://github.com/matplo/heppy) which wraps fastjet and fastjet-contribs in python

#### One-time setup

We recommend to use [pipenv](https://github.com/pypa/pipenv) to manage your python environment:

```
cd JETSCAPE-analysis
pipenv --three
pipenv install pyhepmc_ng pyyaml numpy tqdm ROOT
```

Install `heppy` wherever you desire, and load its modules:

```
cd <my-heppy-location>
git clone git@github.com:matplo/heppy.git
cd heppy
./scripts/setup.sh --buildext --root
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
python analyze_events.py -c ../../config/jetscapeAnalysisConfig.yaml -o /my/outputdir
```

where `/my/outputdir` is the directory containing the generated JETSCAPE events.



### Setup `poetry`

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
