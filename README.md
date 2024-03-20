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
cd jetscape_analysis/generate
python jetscape_events.py -c /home/jetscape-user/JETSCAPE-analysis/config/example.yaml -o /home/jetscape-user/JETSCAPE-analysis-output
```

where
- `-c` specifies a configuration file that should be edited to specify the pt-hat bins and JETSCAPE XML configuration paths,
- `-o` specifies a location where the JETSCAPE output files will be written.

Note that the machinery here only modifies the pt-hat bins and (optionally) other parameter values in the JETSCAPE XML configuration -- but does not allow to change which modules are present -- for that you need to manually edit the user XML file.

That's it! The script will write a separate sub-directory with JETSCAPE events for each pt-hat bin.

## (2) Analyzing events

We provide a simple framework to loop over the generated JETSCAPE output files,
perform physics analysis, and produce a ROOT file.
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

We recommend to use a virtual environment such as [virtualenv]() or [pipenv](https://github.com/pypa/pipenv) to
manage your python environment, e.g.:

```bash
$ cd /home/jetscape-user/JETSCAPE-analysis
$ # Using virtualenv
$ python3 -m venv .venv
$ source .venv/bin/activate
$ pip install pyhepmc pyyaml numpy tqdm
$ # Using pipenv
$ pipenv --three
$ pipenv install pyhepmc pyyaml numpy tqdm
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

The class `jetscape_analysis/analysis/analyze_events_base.py` is a base class to analyze JETSCAPE events and produce an output ROOT file.
To use this, you should write your own class which inherits from `analyze_events_base.py` and implements the following two functions:
- `initialize_user_output_objects()` -- This defines the ROOT histograms or trees that you want to write (called once per JETSCAPE output file)
- `analyze_event(event)` -- This is where you fill your output objects (called once per event). The `event` object and available functions for HepMC or Ascii format can be seen in `jetscape_analysis/analysis/event`.

As an example -- and starting point to copy-paste and build your own class -- see `analyze_events_example.py`.
Simply run the script:

```
cd /home/jetscape-user/JETSCAPE-analysis/jetscape_analysis/analysis
python analyze_events_example.py -c ../../config/example.yaml -i /home/jetscape-user/JETSCAPE-analysis-output -o /my/outputdir
```

where
- `-c` specifies a configuration file that should be edited to specify the pt-hat bins and analysis parameters,
- `-i` specifies  is the directory containing the generated JETSCAPE events,
- `-o` specifies a location where the analysis output will be written.

See `config/example.yaml` for required analysis settings and further details, such as whether to scale and merge the pt-hat bins.

---------------------------------------------------------------------

Please post an [Issue](https://github.com/jdmulligan/JETSCAPE-analysis/issues) if you encounter any problems,
or have suggestions for improvement!

---------------------------------------------------------------------

### Setup package development

This section is only relevant if you want to package up the code yourself. For most other development purposes,
simply working on an editable copy of the repository with `pip install -e .` is sufficient.

To package up the code, we use the `hatchling` build backend. It complies with standards, such that it
interacts well with `pip`, `build`, etc. In order to use a lockfile -- which isn't available via hatch
as of March 2024 -- we use [`pdm`](https://github.com/pdm-project/pdm).

To use `pdm`, you need to set it up once globally. The easiest way to do this is by first installing
[pipx](https://github.com/pypa/pipx), which is broadly available via a package manager for your operating
system. After installing pipx, install `pdm` globally with `pipx install pdm`.

Next, create a virtual environment for this project with:

```bash
$ python3 -m venv .venv
```

and you can enter it with:

```bash
$ source .venv/bin/activate
```

Finally, you can install the package using the lockfile with development packages using:

```bash
$ pdm install -G dev
```

#### Pre-commit checks

To setup checks that run on every commit, run

```bash
$ pipx install pre-commit
$ pre-commit install
```

Now, each commit will be checked on the users' machine.
