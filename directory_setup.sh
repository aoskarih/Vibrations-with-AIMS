#!/bin/bash

# Folder names
# Root folder
root=TCNE

# Subfolders
sub=('relaxations' 'vibrations' 'restart_files')

# State folders are named so that if the name contains "oc" the state needs to force occupations. 
states=('neutral_0' 'ion_oc_0' 'ion_oc_1' 'ion_oc_2' 'ion_oc_3' 'ion_oc_4' 'ion_oc_5' 'ion_oc_6' 'ion_oc_7' 'ion_oc_8' 'ion_oc_9' 'ion_oc_10' 'ion_oc_11' 'ion_oc_12' 'ion_oc_13')

# Script names
run_vib=run_vib.sh
run_rel=run_relax.sh
get_vib=get_vibrations.py
get_occ=get_vibrations_occ.py

# Chosen delta
vib_run_dir=delta_0.0025

cur=$PWD

mkdir $root

cp FCI.py $root
cp $run_vib $root
cp $run_rel $root
cp $get_vib $root
cp $get_occ $root

cd $root

 
tmp=${sub[0]}
mkdir $tmp
for j in "${states[@]}"; do
    mkdir $tmp/$j
    cp ${run_rel} $tmp/$j/
done

tmp=${sub[1]}
mkdir $tmp
for j in "${states[@]}"; do
    mkdir $tmp/$j
    re='(oc)'
    if [[ $j =~ $re ]]; then
        sed "s+get_vibrations.py+get_vibrations_occ.py+" $run_vib > tmp_file
        mv tmp_file $tmp/$j/$run_vib
        sed "s+path_to_restart_files+$cur/$root/${sub[2]}/$j/$vib_run_dir/+" $get_occ > tmp_file
        mv tmp_file $tmp/$j/$get_occ
    else
        cp $get_vib $tmp/$j/
        cp $run_vib $tmp/$j/
    fi
done

tmp=${sub[2]}
mkdir $tmp
for j in "${states[@]}"; do
    mkdir $tmp/$j
    cp $get_vib $tmp/$j/
    cp $run_vib $tmp/$j/
done
 



