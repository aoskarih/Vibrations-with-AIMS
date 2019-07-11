#!/bin/bash

# Folder names
root=Benzene
sub=('relaxations' 'vibrations' 'restart_files')
states=('neutral' 'ion_0' 'ion_1' 'ion_2' 'ion_3' 'ion_4')



mkdir $root
cd $root

for i in "${sub[@]}"; do
    mkdir $i
    cd $i
    for j in "${states[@]}"; do
        mkdir $j
    done
    cd ..
done
 



