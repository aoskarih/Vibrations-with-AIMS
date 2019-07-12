
# Vibrational spectrum with FHI-AIMS


List of contents
----------------







Making spectrum for molecule
----------------------------

This is a step-by-step guide to getting vibrational spectrum for your molecule.

### Steps

1. [Setting up directory](#setting-up-directory)
2. [Running relaxation calculations](#running-relaxation-calculations)
3. [Running ground state vibrational calculations](#running-ground-state-vibrational-calculations)
4. [Running exited state vibrational calculations](#running-exited-state-vibrational-calculations)
5. [Post-processing data from calculations](#post-processing-data-from-calculations)
6. [Plotting results](#plotting-results)
7. [Summary](#summary)

### Setting up directory

Some of the scripts require specific folder structure and folder names, so if you make changes please be sure to check that scripts still work. There should be place in scripts to easily edit the place where it looks for files.

Benzene will be used as an example throughout this guide and all the default setting are for benzene. This means that we will have neutral ground state from which all the transitions are going to start and we will take five ion states where transitions are going to end. Ion states are going to differ in their electronic occupations.

#### Folders

The script `directory_setup.sh` will make the folders and copy files automatically, but first you should check that scripts `run_relax.sh` and `run_vib.sh` really run AIMS in your system and that list of states and folder names in `directory_setup.sh` match your needs.

First a root folder should be made and three sub folder inside it. The subfolders are for relaxation calculations, vibrational calculations and restart files. Restart files are needed when forcing occupations. Each of these three is also divided to subfolders for each state. Now in this example we have six of them. Now the directory should look something like below.
```
Benzene/
    relaxations/
        neutral/
        ion_0/
        ion_1/
        ion_2/
        ion_3/
        ion_4/
    vibrations/
        neutral/
        ...
    restart_files/
        neutral/
        ...
```

#### Files

To run AIMS we need different script depending on the folder. For the "*relaxations*" folder `run_relax.sh` is used. For states in folders "*vibrations*" and "*restart_files*", `run_vib.sh` is used. `run_vib.sh` manages folders and runs `get_vibrations_occ.py` if state must force occupations and `get_vibrations.py` if not. For more details about `get_vibrations.py` see AIMS Manual section "_Calculation of vibrational and phonon frequencies_". `get_vibrations_occ.py` is slightly modified version of the same script and only change as of writing this is at line 350.

### Running relaxation calculations

To run any calculation we of course need `control.in` and `geometry.in` files. Sample files for benzene are provided. 

### Running ground state vibrational calculations


delta in 2 places


### Running exited state vibrational calculations


### Post-processing data from calculations


### Plotting results


### Summary


