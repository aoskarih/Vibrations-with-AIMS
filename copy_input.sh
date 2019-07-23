#!/bin/bash

initial=neutral_0

for d in relaxations/*/ ; do
    cp relaxations/$initial/control.in restart_files/$d
    cp relaxations/$d/geometry.in.next_step restart_files/$d/geometry.in
    cp relaxations/$d/control.in vibrations/$d
    cp relaxations/$d/geometry.in.next_step vibrations/$d/geometry.in
done
