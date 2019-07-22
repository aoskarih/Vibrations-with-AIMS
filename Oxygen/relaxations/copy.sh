#!/bin/bash

for d in */ ; do
    re="(neut_0)"
    if [[ $d == $re ]]; then
        echo "no $d"
    else
        echo $d
        cp $d/control.in ../vibrations/$d
        cp $d/geometry.in.next_step ../vibrations/$d/geometry.in
    fi
done

