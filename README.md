
# Vibrational spectrum with FHI-AIMS

### Steps

1. [Setting up directory](#setting-up-directory)
2. [Running relaxation calculations](#running-relaxation-calculations)
3. [Generating restart files](#generating-restart-files)
4. [Running vibrational calculations](#running-vibrational-calculations)
5. [Calculating transition intensities](#calculating-transition-intensities)
6. [Plotting results](#plotting-results)

A. [Summary](#summary)\
B. [List of scripts](#list-of-scripts)

### Setting up directory

Some of the scripts require specific folder structure and folder names, so if you make changes please be sure to check that scripts still work. There should be place in scripts to easily edit the place where it looks for files.

Benzene will be used as an example throughout this guide and all the default setting are for benzene. This means that we will have neutral ground state from which all the transitions are going to start and we will take five ion states where transitions are going to end. Ion states are going to differ in their electronic occupations.

The script `directory_setup.sh` will make the folders and copy files automatically, but first you should check that scripts `run_relax.sh` and `run_vib.sh` really run AIMS in your system and that list of states and folder names in `directory_setup.sh` match your needs. In the script `run_vib.sh` there are two variables *aims_bin* for the location of AIMS binary and *suffix* for the suffix used when running AIMS. The suffix should be put into quotes e.g. `suffix='mpirun -np 4'`.

Now is also the best time to choose the delta that is used as a displacement in the vibrational calculations. The default is 0.0025 Å and it can be changed with script. For example to set *delta=0.001* `./set_delta.sh 0.001`.

#### Folders

First a root folder should be made and three sub folder inside it. The subfolders are for relaxation calculations, vibrational calculations and restart files. Restart files are needed when forcing occupations. Each of these three is also divided to subfolders for each state. Now in this example we have six of them. Now the directory should look something like below.
```
Benzene/
    relaxations/
        neutral_0/
        ion_oc_0/
        ion_oc_1/
        ion_oc_2/
        ion_oc_3/
        ion_oc_4/
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

#### States without `force_occupation_projector`

To run any calculation we need `control.in` and `geometry.in` files. Sample files for benzene are provided. 

When making control.in file it is important to include the line `restart        restart` because the restart file is needed later when forcing the occupations. For more details about forcing the occupations see tag "*force_occupation_projecto*" in AIMS Manual. Name of the restart file should be "restart" and if you wan't to change it, then the edited line (350) in `get_vibrations_occ.py` should also be changed. It is also important to remember that you can only change charge and add the occupation forcing in `control.in` file after running the first relaxation. 

When the `control.in` and `geometry.in` files are ready for the neutral state (or any other initial state) with no occupations, you should copy files to the appropriate folder, in our example case `/Benzene/relaxations/neutral_0/` and then run AIMS. This should be done for all the states that don't require forcing occupations. With benzenes case `ion_0` is the only other.

#### States with `force_occupation_projector`

After AIMS has done its work there should be `restart` file in the folder. Copy that file and the `geometry.in` file to folders of all the states that force occupations and copy `control.in` file to the same folders. Then add the "*force_occupation_projector*" lines for different states. Also remember to change charge for ions. With our example following lines are added
```
ion_0/control.in: force_occupation_projector      21 1 1.0 20 21
ion_1/control.in: force_occupation_projector      20 1 1.0 20 21
ion_2/control.in: force_occupation_projector      19 1 1.0 18 19
ion_3/control.in: force_occupation_projector      18 1 1.0 18 19
ion_4/control.in: force_occupation_projector      17 1 1.0 17 17
```
Then run AIMS for all the states. After AIMS has finished there should be `geometry.in.next_step` file in every folder. Those are needed in the next step.

There might be problems when relaxing the geometries. Try running `grep "Maximum force" */aims.out` in "*relaxations*" folder to make sure all the geometries are converged. Trying to relax all the geometries is the most inconsistent step, especially if there's a lot of states close to each other in energy. Usually after enough playing around with the last two numbers of `force_occupation_projector` all the geoetries converge.

If you have to force occupations in degenerate electron levels, be aware that AIMS usually chooses one of the states and you can't force hole to the other without running a risk that calculation doesn't converge. This means that for example in benzenes case *ion_2* and *ion_3* states will be identical as well as *ion_0* and *ion_1*.

### Generating restart files

Because calculations, where occupations are forced, need restart files and for restart file to work `geometry.in` files must be identical, we need 6N restart files per state (N is number of atoms). These files can be created by running "dummy" vibrational calculations with the control file of the ground state and the geometry file of the exited state.

To run these calculations, first copy `geometry.in.next_step` files from "*relaxations*" to matching folders in "*restart_files*" and rename them `geometry.in`. Then copy the `control.in` file that was used for the initial state to all of the folders in "*restart_files*". You can also run the script `copy_input.sh` which will copy the files. It should be ran in root directory ("*Benzene*" in our example).

Restart files are needed only for the states that force occupations, so for our example we don't need to touch *neutral_0* state at all. But the other states are ready and next `run_vib.sh` should be ran in those folders. It will make new folder "*delta_0.0025*" and vibrational calculations will be done there. Because `control.in` and `geometry.in` files don't match in these calculations, results will be wrong. But restart files will be generated for all the steps of the calculation and these will be used in the next step.

After the script has finished there should be list of folders under "*delta_0.0025*" with names something like "*run.i_atom_1.i_coord_1.displ_0.0025*" and inside everyone of these should be a restart file.

The 0.0025 in "*delta_0.0025*" is the displacement used when doing the vibrational calculations. If you want to change the value of delta after generating the restart files, then `get_vibrations_occ.py` and `run_vib.sh` scripts must be edited in every subfolder as well as `FCI.py` in root folder. Remember that restart files and final calculations must be done with the same delta.

### Running vibrational calculations

Now that we have all the restart files and proper geometries we can begin the actual vibrational calculations. First we need to copy relaxed geometries and `control.in` files from "*relaxations*" to "*vibrations*". This time we want the `control.in` files specific to state and again remember to rename `geometry.in.next_step` to just `geometry.in`. If you ran the `copy_input.sh` previously, then all the files shoudl be copied already. When these files are in place, everything should be ready and you only need to run the script `run_vib.sh` in every directory under "*vibrations*". 

After aims has finished there will be output files in all of the `vibrations/state/delta_0.0025/` directories. What we are interested in are the `run.xyz` and `vib_post_0.0025.out` files. You should check that values in `run.xyz` are in float format. Some times they appear as complex numbers, example: `0.12345+0.00000j`. That can be usually fixed by running the script again: `python get_vibrations.py run 1 >& vib_post_0.0025.out`. 

In addition to normal modes, there's values for force constants, reduced masses and frequencies in `run.xyz`. Sometimes the value for force constant or reduced mass is replaced by \*\*\*\*\*. Post-processing script can handle the situations when only one of the values is missing, but when both of them are missing it causes problems. 

It should also be checked that there isn't any unstable frequencies. When calculating the transitions most of the modes are dropped out, so if there's one or two unstable frequensices it shouldn't affect the results that much, but if one of the relevant modes is missing, it could have noticeable effect.

### Calculating transition intensities

When the output files of the vibrational calculations are ready, first should be checked that the folder and filenames are correct at the top of the script `FCI.py`. Then the script should be ran in root directory which for our example is "*Benzene*". The script will make a `intensity.dat` file which will have the spectrum we are after. You can change the filename at the bottom of `FCI.py`. When running the script it is recommended to direct the output to a file, for example `python FCI.py >& fci.out`. 

There's few things that can be modified in `FCI.py`. At the top of the file there's variables `S_lim`, `m_lim` and `n_lim`. These control the number of transitions that is calculated. `S_lim` drops out vibrational modes that aren't relevant. `m_lim` is the ceiling for the sum of vibrational quantum numbers of the initial state and `n_lim` is the corresponding value for the final state. There shouln't be any problems with the runtime of the script as long as `m_lim` and `n_lim` are kept reasonably small (defaults: *m=2*, *n=3*). If script seems to take too long to finish even with small *m* and *n* check the output of the script to see if there's a large number (>10) of relevant modes. If so, then make `S_lim` larger.

The script can also be given a file which sets the positions of 0-0 peaks. If no such file is provided the script will read total energies from the relaxation calculations and calculate the positions for the peaks from those values. The file should have the energy of initial state in the first line and then the energies of the final states on their own lines in the same order that they are calculated in the script. Example:
```
0.000
9.055
9.498
12.304
12.318
12.458
```
Now the first value is zero, but only thing that matters is the differences in the values. There can be more values in the file than what is needed as long as the needed values are first. In general it is best to use file for the 0-0 positions since the values from relaxations aren't that good and the script fails if any of the relaxations has failed.

The output file is a list of energies and corresponding intensities for different transitions. There's also identifier which tells the final state associated with the transition.

### Plotting results

There's a script `plot_intensity.py` which can plot the data that `FCI.py` outputs. It can also calculate the shift from 0-0 peak caused by vibrations. In the bottom of the `plot_intensities.py` there's place where functions can be called to make plots.

`plot_data(...)` plots the delta peaks and does the gaussian broadening to them. `plot_peak_shift()` calculates and draws the shift in peak intensity caused by vibrations. And here's the results for our benzene example, first `plot_data()`
![Energy spectra of benzene](https://raw.githubusercontent.com/aoskarih/Vibrations-with-AIMS/master/example_files/benzene_0-4_GW.png)
And then `plot_peak_shift()`
![Energy spectra of benzene](https://raw.githubusercontent.com/aoskarih/Vibrations-with-AIMS/master/example_files/peak_shift.png)


### Summary

This is a summary of the steps 1-5 assuming all default settings and no problems.

1. *Setup*\
By editing `directory_setup.sh`, change the name of the root directory to your liking and then define all the states that you want in your calculations. Then run the script.

2. *Relaxation*\
Make `control.in` files for all of your states and make one `geometry.in` to begin with. Copy these to their respective folders and run the relaxations without occupations. When AIMS is ready, copy `restart` file and run the relaxations with occupations.

3. *Restart files*\
First run `copy_input.sh` in root directory. Then run the calculations in "*restart_files*".

4. *Vibrational calculations*\
Run the calculations in "*vibrations*" directory.

5. *Transitions*\
Check that settings are correct in `FCI.py` and then run it.

### List of scripts

+ `copy_input.sh`\
Script handles the copying required in the steps 3 and 4.

+ `directory_setup.sh`\
Script sets up folders and copies some necessary files to them. It also modifies the `get_vibrations_occ.py` scripts to copy right restart files when ran. Variable `root` and list `states` should be checked and modified before running the script.

+ `FCI.py`\
Script calculates the transitions and their intensities between electronic states that are defined in the list `folders` in the beginning of the file. Transitions start always from the same electronic state which is the first state in the list. 

+ `get_vibrations.py`\
Script for calculating vibrational modes and frequencies. For more info see section "*Calculation of vibrational and phonon frequencies*" in AIMS manual.

+ `get_vibrations_occ.py`\
Slightly modified version of `get_vibrations.py`. Copies the restart files necessary for calculations with occupations.

+ `plot_intensity.py`\
Script for plotting the output of `FCI.py`. 

+ `run_relax.sh`\
Runs AIMS and directs output to a file. Make sure that script really finds AIMS.

+ `run_vib.sh`\
Sets up a directory and runs vibrational calculations. Make sure that script really finds AIMS.

+ `set_delta.sh`\
If you want to change the delta used in the vibrational calculations, then run this script with the new delta as a parameter before first step. 

