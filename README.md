
# Vibrational spectrum with FHI-AIMS


## Making vibrational spectrum

This is a step-by-step guide to getting vibrational spectrum for your molecule with FHI-AIMS.

### Steps

1. [Setting up directory](#setting-up-directory)
2. [Running relaxation calculations](#running-relaxation-calculations)
3. [Generating restart files](#generating-restart-files)
4. [Running vibrational calculations](#running-vibrational-calculations)
5. [Calculating transition intensities](#calculating-transition-intensities)
6. [Plotting results](#plotting-results)
7. [Summary](#summary)

### Setting up directory

Some of the scripts require specific folder structure and folder names, so if you make changes please be sure to check that scripts still work. There should be place in scripts to easily edit the place where it looks for files.

Benzene will be used as an example throughout this guide and all the default setting are for benzene. This means that we will have neutral ground state from which all the transitions are going to start and we will take five ion states where transitions are going to end. Ion states are going to differ in their electronic occupations.

The script `directory_setup.sh` will make the folders and copy files automatically, but first you should check that scripts `run_relax.sh` and `run_vib.sh` really run AIMS in your system and that list of states and folder names in `directory_setup.sh` match your needs. 

#### Folders

First a root folder should be made and three sub folder inside it. The subfolders are for relaxation calculations, vibrational calculations and restart files. Restart files are needed when forcing occupations. Each of these three is also divided to subfolders for each state. Now in this example we have six of them. Now the directory should look something like below.
```
Benzene/
    relaxations/
        neutral_0/
        ion_0/
        ion_1/
        ion_2/
        ion_3/
        ion_4/
    vibrations/
        neutral_0/
        ...
    restart_files/
        neutral_0/
        ...
```

#### Files

To run AIMS we need different script depending on the folder. For the "*relaxations*" folder `run_relax.sh` is used. For states in folders "*vibrations*" and "*restart_files*", `run_vib.sh` is used. `run_vib.sh` manages folders and runs `get_vibrations_occ.py` if state must force occupations and `get_vibrations.py` if not. For more details about `get_vibrations.py` see AIMS Manual section "*Calculation of vibrational and phonon frequencies*". `get_vibrations_occ.py` is slightly modified version of the same script and the only change as of writing this is at line 350.

### Running relaxation calculations

To run any calculation we need `control.in` and `geometry.in` files. Sample files for benzene are provided. 

When making control.in file it is important to include the line `restart        restart` because the restart file is needed later when forcing the occupations. For more details about forcing the occupations see tag "*force_occupation_projecto*" in AIMS Manual. Name of the restart file should be "restart" and if you wan't to change it, then the edited line (350) in `get_vibrations_occ.py` should also be changed. It is also important to remember that you can only change charge and add the occupation forcing in `control.in` file after running the first relaxation. 

When the `control.in` and `geometry.in` files are ready for the neutral state (or any other initial state) with no occupations, you should copy files to the appropriate folder, in our example case `/Benzene/relaxations/neutral_0/` and then run AIMS. This should be done for all the states that don't require forcing occupations. With benzenes case `ion_0` is the only other.

After AIMS has done its work there should be `restart` file in the folder. Copy that file and the `geometry.in` file to folders of all the states that force occupations and copy `control.in` file to the same folders. Then add the "*force_occupation_projector*" lines for different states. Also remember to change charge for ions. With our example following lines are added
```
ion_1/control.in: force_occupation_projector      20 2 0.0 20 21
ion_2/control.in: force_occupation_projector      19 2 0.0 18 19
ion_3/control.in: force_occupation_projector      18 2 0.0 18 19
ion_4/control.in: force_occupation_projector      17 2 0.0 16 18
```
Then run AIMS for all the states. After AIMS has finished there should be `geometry.in.next_step` file in every folder. Those are needed in the next step.

There might be problems when relaxing the geometries. Try running `grep "Maximum force" */aims.out` in "*relaxations*" folder to make sure all the geometries are converged.

If you have to force occupations in degenerate electron levels, be aware that AIMS usually chooses one of the states and you can't force hole to the other without running a risk that calculation doesn't converge. This means that for example in benzenes case *ion_2* and *ion_3* states will be identical as well as *ion_0* and *ion_1*.

### Generating restart files

Because calculations, where occupations are forced, need restart files and for restart file to work `geometry.in` files must be identical, we need 6N restart files per state (N is number of atoms). To make these files, first copy `geometry.in.next_step` files from "*relaxations*" to matching folders in "*restart_files*" and rename them `geometry.in`. Then copy the `control.in` file that was used for the initial state to all of the folders in "*restart_files*". 

Restart files are needed only for the states that force occupations, so for our example we don't need to touch *neutral_0* or *ion_0* states at all. But the other states are ready and next `run_vib.sh` should be ran in those folders. It will make new folder "*delta_0.0025*" and vibrational calculations will be done there. Because `control.in` and `geometry.in` files don't match in these calculations, results will be wrong. But restart files will be generated for all the steps of the calculation and these will be used in the next step.

After the script has finished there should be list of folders under "*delta_0.0025*" with names something like "*run.i_atom_1.i_coord_1.displ_0.0025*" and inside everyone of these should be a restart file.

The 0.0025 in "*delta_0.0025*" is the displacement used when doing the vibrational calculations. If you want to change the value of delta, it has to be changed in three places: `run_vib.sh`, `directory_setup.sh` and `FCI.py`. Script `set_delta.sh` does this automatically and new delta is the first parameter, for example `./set_delta.sh 0.001`. Remember that restart files and final calculations must be done with the same delta.

### Running vibrational calculations

Now that we have all the restart files and proper geometries we can begin the actual vibrational calculations. First we need to copy relaxed geometries and `control.in` files from "*relaxations*" to "*vibrations*". Again remember to rename `geometry.in.next_step` to just `geometry.in`. When these files are in place, everything should be ready and you only need to run the script `run_vib.sh` in every directory under "*vibrations*". 

After aims has finished there will be output files in all of the `vibrations/state/delta_0.0025/` directories. What we are interested in are the `run.xyz` and `vib_post_0.0025.out` files. You should check that values in `run.xyz` are in float format. Some times they appear as complex numbers, example: `0.12345+0.00000j`. That can be usually fixed by running the script again: `python get_vibrations.py run 1 >& vib_post_0.0025.out`. 

In addition to normal modes, there's values for force constants, reduced masses and frequencies in `run.xyz`. Sometimes the value for force constant or reduced mass is replaced by \*\*\*\*\*. Post-processing script can handle the situations when only one of the values is missing, but when both of them are missing it causes problems.

### Calculating transition intensities

When the output files of the vibrational calculations are ready, the script `FCI.py` should be ran in root directory which for our example is "*Benzene/*". The script will make a `intensity.dat` file which will have the spectrum we are after. You can change the filename at the bottom of `FCI.py`. When running the script it is recommended to direct the output to a file, for example `python FCI.py >& fci.out`. 



### Plotting results


### Summary



