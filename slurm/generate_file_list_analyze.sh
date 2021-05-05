#! /usr/bin/bash

for d in */ ; do
    cd $d
    find "$PWD" -name "JetscapeHadron*" > files.txt
    cd ..
done
