#! /bin/bash

PREFIX="/home/jetscape-user/JETSCAPE-analysis"
OUTPUT_DIR="test/20200219"
REFERENCE_DIR="test/20200218"

echo ""
echo "Running Tests..."
echo ""

cd $PREFIX/jetscape_analysis/generate
python jetscape_events.py -c $PREFIX/config/jetscapeTestConfig.yaml -o $PREFIX/$OUTPUT_DIR

echo ""
echo "Comparing tests in $OUTPUT_DIR to Reference in $REFERENCE_DIR ..."
echo ""

N=0
N_PASSED_HEPMC=0
N_PASSED_ASCII=0
cd $PREFIX/$OUTPUT_DIR
for dir in */ ; do
  N=$((N+1))
  
  DIFF_HEPMC=$(diff $PREFIX/$OUTPUT_DIR/$dir/test_out.hepmc $PREFIX/$REFERENCE_DIR/${dir}test_out.hepmc)
  if [ $? -ne 0 ]
  then
    echo "Error: Check whether you have used the same YAML config for the Test and the Reference"
    exit 1
  fi
  
  DIFF_ASCII=$(diff $PREFIX/$OUTPUT_DIR/$dir/test_out.dat $PREFIX/$REFERENCE_DIR/${dir}test_out.dat)
  if [ $? -ne 0 ]
  then
    echo "Error: Check whether you have used the same YAML config for the Test and the Reference"
    exit 1
  fi

  if [ "${DIFF_HEPMC}" == "" ]
  then
    N_PASSED_HEPMC=$((${N_PASSED_HEPMC}+1))
  else
    echo "Test $dir failed for HepMC"
  fi

  if [ "${DIFF_ASCII}" == "" ]
  then
    N_PASSED_ASCII=$((${N_PASSED_ASCII}+1))
  else
    echo "Test $dir failed for Ascii"
  fi
    
done

N_FAILED_HEPMC=$(($N-$N_PASSED_HEPMC))
N_FAILED_ASCII=$(($N-$N_PASSED_ASCII))
if [[ $N_FAILED_HEPMC -eq 0 && $N_FAILED_ASCII -eq 0 ]]
then
  echo "All $N tests passed! :)"
else
  echo ""
  echo "Tests FAILED :("
  echo "$N_FAILED_HEPMC/$N tests FAILED for HepMC"
  echo "$N_FAILED_ASCII/$N tests FAILED for Ascii"
  exit 1
fi
echo ""


