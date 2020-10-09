#! /usr/bin/bash

for d in */ ; do
    cd $d
    find "$PWD" -name "JetscapeHadron*.out" > files.txt
    cd ..
done
