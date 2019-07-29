#!/bin/bash

for d in */ ; do
    re="(neutra)"
    if [[ $d =~ $re ]] ; then
        echo no $d
    else
        cd $d
        sbatch $1
        cd ..
    fi
done
