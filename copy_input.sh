#!/bin/bash

initial=neutral_0
cd relaxations
for d in */ ; do
    d=${d%*/}
    echo $d
    cp $initial/control.in ../restart_files/$d
    cp $d/geometry.in.next_step ../restart_files/$d/geometry.in
    cp $d/control.in ../vibrations/$d
    cp $d/geometry.in.next_step ../vibrations/$d/geometry.in
done
