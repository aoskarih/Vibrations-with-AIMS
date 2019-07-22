#!/bin/bash

for d in */ ; do
    re="(_0)"
    if [[ $d =~ $re ]]; then
        echo "no $d"
    else
        echo $d
        cp ion_oc_0/run_vib.sh $d
    fi
done

