
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


### Setting up directory

Some of the scripts require specific folder structure and folder names, so if you make changes please be sure to check that scripts still work. There should be place in scripts to easily edit the place where it looks for files.

Benzene will be used as an example throughout this guide. This means that we will have neutral ground state from which all the transitions are going to start and we will take five ion states where transitions are going to end. Ion states are going to differ in their electronic occupations.

#### Folders

The script `directory_setup.sh` will do this step automatically. You should also check that scripts `run_relax.sh`, `get_vibrations.py` and `get_vibrations_occ.py` really run AIMS in your system.

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

To run AIMS we need different script depending on the folder. For the relaxations `run_relax.sh` is used. This is optional and it is used only for conveniency. For states in folder "*vibrations*" `get_vibrations_occ.py` is used if state must force occupation and `get_vibrations.py` if not. In folder "*restart_files*" only `get_vibrations.py` is used.



### Running relaxation calculations


### Running ground state vibrational calculations


### Running exited state vibrational calculations


### Post-processing data from calculations


### Plotting results




