## Overview

This repository contains basic tools to run and analyze Jetscape events in HepMC3 format.

The script `doJetAnalysis.py` is the top-level script, which allows you to:
1. Generate Jetscape events, including automated machinery to launch a set of pt-hat bins
2. Analyze Jetscape events, producing an output ROOT file
3. Aggregate the results from the set of pt-hat bins, and plot the analysis results

The analysis machinery (step 2) consists of a run macro `runJetscapeAnalysis.C` which steers
an analysis task `JetscapeAnalysis`.

## Pre-requisites
These tools exist independently of the Jetscape framework, however
it is recommended to run the script from inside the Jetscape docker container,
in order to automatically satisfy pre-requisites (HepMC, ROOT).
Assuming you have a Jetscape docker according to the Jetscape instructions 
(with a shared folder located at `~/jetscape-docker`, containing the Jetscape repository at `~/jetscape-docker/JETSCAPE`), 
you should do (from outside the docker container):

```
cd ~/jetscape-docker/
git clone git@github.com:jdmulligan/JETSCAPE-analysis.git
```

And then use the same workflow as with Jetscape: Build and run the analysis machinery inside
the docker container, but access the source and output files from outside the container.

## Building
In order to perform the analysis (step 2), you must build using cmake:
```
mkdir build
cd build
cmake ..
make
```

This will produce an executable `runJetscapeAnalysis`, which can be run via `doJetAnalysis.py`,
or it can be run directly on a hepmc output file: `./runJetscapeAnalysis myOutput.hepmc`.
