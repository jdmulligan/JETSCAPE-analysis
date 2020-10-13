#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

# Command line arguments
if [ "$1" != "" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  JOB_ID=$2
  echo "Job ID: $JOB_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  TASK_ID=$3
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

# Define output path from relevant sub-path of input file
# Note: suffix depends on file structure of input file -- need to edit appropriately
OUTPUT_SUFFIX1=$(echo $INPUT_FILE | cut -d/ -f5-5)
FILENAME=$(echo $INPUT_FILE | cut -d/ -f6-6)
OUTPUT_SUFFIX2=$(echo $FILENAME | cut -d. -f 1)
OUTPUT_DIR="/rstorage/jetscape/JETSCAPE-AA-events/skim/$JOB_ID/$OUTPUT_SUFFIX1/$OUTPUT_SUFFIX2"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list

# Run python script via pipenv
cd /software/users/james/jetscape-docker/JETSCAPE-analysis
pipenv run python jetscape_analysis/analysis/reader/skim_ascii.py -i $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /rstorage/jetscape/JETSCAPE-AA-events/skim/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/jetscape/JETSCAPE-AA-events/skim/${JOB_ID}/
