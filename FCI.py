# Arttu Hyvonen
# 7/2019

import scipy.special as special
from scipy.misc import factorial
import numpy as np
import matplotlib.pyplot as plt
from  heapq import nlargest

# constants
hbar = 6.5821e-16	# eV s
c = 2.9979e8		# m/s
pi = 3.1416
u = 1.6605e-27		# kg
kb = 8.6173e-5		# ev/K


# Molecule specific settigs
##################################

# Limiting number of transitions
m_lim = 2
n_lim = 3
S_lim = 1e-4


delta = 0.0025

# data folders
# Sub folder where output files are
sub_fol = "delta_"+str(delta)+"/"          

# List of folders for different electronic states
# Transitions are calculated from the first state to the others
folders = ["neutral_0/",                    
            "ion_0_oc/",                       
            "ion_1_oc/",  
            "ion_2_oc/", 
            "ion_3_oc/", 
            "ion_4_oc/"]

# file names
xyz = "run.xyz"                 # xyz file which contains normal modes
vib = "vib_post_0.0025.out"     # output file which contains frequencies in a list
mas = "masses.run_0.0025.dat"   # mass file, currently data from here not in use

output_file = "intensity.dat"   # Output filename
energy_file = ""                # Filename for the 0-0 peak energies. If left empty script will use energies from relaxations.

##################################


# returns positions of the atoms from the .xyz file
# Angstroms
def get_positions(fold):
    positions = []
    for folder in fold:
        f = open(folder+xyz, "r")
        n = int(next(f).split()[0])
        next(f)
        pos = np.zeros(3*n)
        for i in range(n):
            l = next(f).split()
            for j in range(3):
                pos[3*i+j] = float(l[j+1])
        positions.append(pos)
        f.close()
    return positions

# returns normal modes from the .xyz file
def get_normal_modes(fold):
    mode_lists = []
    for folder in fold:
        f = open(folder+xyz, "r")
        n = int(next(f).split()[0])
        f.seek(0)
        modes = []
        for _ in range(3*n):
            mode = np.zeros(3*n)
            next(f)
            next(f)
            for i in range(n):
                l = next(f).split()
                for j in range(3):
                    mode[3*i+j] = float(l[j+4])
            modes.append(mode)
        mode_lists.append(modes)
        f.close()
    return mode_lists

# returns frequencies from the .out file
# 1/s
def get_frequencies(fold):
    frequency_lists = []
    for folder in fold:
        f = open(folder+vib, "r")
        n = 0
        for l in f:
            if "Number of atoms" in l:
                n = int(l.split()[4])
                continue
            if "Mode number" in l:
                break
        freq_list = np.zeros(3*n)
        for i in range(3*n):
            freq = float(next(f).split()[1])
            freq *= 100*c*2*pi			# 1/cm -> 1/s
            freq_list[i] = freq
        frequency_lists.append(freq_list)
    return frequency_lists

# returns transfromation matrix from cartesian coordinates to coordinates in given basis 
# basis vectors are assumed to be in cartesian coordinates
def get_transformation_matrix(basis_vectors):
    n = len(basis_vectors)
    T = np.zeros((n,n))
    for i in range(n):
        v = basis_vectors[i]
        l = np.linalg.norm(v)
        v = v/l
        for j in range(n):
            T[j][i] = v[j]
    return np.linalg.inv(T)

# returns reduced masses for the different vibrational modes from .xyz
# kg
def get_reduced_masses(fold):
    list_list = []
    for folder in fold:
        f = open(folder+xyz, "r")
        mass_list = []
        for l in f:
            if "stable frequency" in l:
		if l.split()[13] == "*****" and l.split()[18] == "*****":
		    mass = 0
                elif l.split()[13] == "*****":
		    mass = float(l.split()[18])/(float(l.split()[3])*100*c*2*pi)**2
	        else:
		    mass = float(l.split()[13])
		mass_list.append(mass*u)	# u -> kg
        list_list.append(mass_list)
    return list_list

# returns force constants (k) for the different vibrational modes from .xyz
# kg/s^2
def get_force_constants(fold):
    list_list = []
    for folder in fold:
        f = open(folder+xyz, "r")
        k_list = []
        for l in f:
            if "stable frequency" in l:
	        if l.split()[13] == "*****" and l.split()[18] == "*****":
		    k = 0
		elif l.split()[18] == "*****":
		    k = float(l.split()[13])*u*(200*c*pi*float(l.split()[3]))**2
		else:
		    k = float(l.split()[18])
		k_list.append(k*100)		# mDyne/A -> kg/s^2
        list_list.append(k_list)
    return list_list

# returns masses from the masses.* file
# kg
def get_masses(fold):
    f = open(fold[0]+mas, "r")
    masses = []
    for l in f:
        m = float(l.split()[0])
        masses.append(m*u)			# u -> kg
    return masses

# returns total energies (no vibration energy) of neutral and cation molecule
# eV
def get_total_energies():
    f = open("total_energies.dat", "r")
    e = []
    for l in f:
	e.append(float(l.split()[0]))
    return e

# returns value of Huang-Rhys parameter
def get_huang_rhys(d, mu, f):
    S = d*d*mu*f/(2*hbar)
    # unit conversion
    S *= 0.062415				# A^2 * kg / (eV * s^2) = 0.0624150...
    return S

# returns value of Franck-Condon integral squared
def FCI2(S, m, n):
    I = np.exp(-S)*S**(n-m)*factorial(m)/factorial(n)*laguerre(S, m, n-m)**2
    return I

# returns value of laguerre polynomial
def laguerre(x, n, k):
    if k > -1:
        L = special.genlaguerre(n, k)
        return L(x)
    else:
        # scipy doesn't handle negative k so recursive equation is used
        if n == 0:
            return 1
        else:
            return laguerre(x, n, k+1) - laguerre(x, n-1, k+1)

# returns relative intensity of vibrational transition
# by passing list relative, which contains modes for which S is above some limit you can make calculation much faster
# example relevant = [9, 17, 21, 26, 31]
def intensity(d, freq, mu, m, n, T, S_limit, relevant=[]):
    I = 1
    if relevant:
        # Do the calculation for all relevant modes
        for i in relevant:
            # Get S
	    S = get_huang_rhys(d[i], mu[i], freq[i])
	    
            # Calculate Franck-Condon integral
            FCI = FCI2(S, m[i], n[i])
            
            # Intensity is integrals for all modes multiplied together
	    I *= FCI
	    
            # Temperature term
            if T != 0:
		I *= np.exp(-hbar*m[i]*freq[i]/(kb*T))

        if I > 1e-3:
            print("I: " + str(I))
	return I
    
    # If no relevant list get S for all modes and continue calculation with relevant modes
    for i in range(len(d)):
	S = get_huang_rhys(d[i], mu[i], freq[i])
        if S < S_limit: 
            #print("%2.2i  %10.6f  %10s %5.1i %5.1i  ignored" % (i, S, "", m[i], n[i]))
            continue
        FCI = FCI2(S, m[i], n[i])
        I *= FCI
	if T != 0:
	    I *= np.exp(-hbar*m[i]*freq[i]/(kb*T))
        #print("%2.2i  %10.6f  %10.6f %5.1i %5.1i" % (i, S, FCI, m[i], n[i]))
    print("I: " + str(I) + "\n")
    return I

# returns energy of some state with m = [...]
# E0 is the zero-point vibrational energy
# eV
def get_vib_energy(freq, m, E0):
    E = E0
    for i, f in enumerate(freq):
        E += f*hbar*m[i]
    return E

# returns list of different integer combinations of length l 
# for which the sum of integers in the combination is less than m.
# Note that the length of the list grows rapidly: length = Binomial(l+m, l)
# example:
# m = 2, l = 4
# returns: [[0, 0, 0, 0], [1, 0, 0, 0]. [2, 0, 0, 0], [1, 1, 0, 0], ..., [0, 0, 0, 2]]
def combinations(m, l):
    arr = []
    ar = [0 for _ in range(l)]
    i = 0
    while i < l:
        if ar not in arr:
            arr.append(ar[:])
        ar[i] += 1
        s = sum(ar)
        if s < m + 1:
            i = 0
        else:
            i += 1
            for j in range(i):
                ar[j] = 0
    return arr


# Writes intensities to file specified with "filename"
def write_intensities(d, freq, mu, E0, T, S_limit, c, filename):
    N = len(d)
    m = np.zeros(N)
    n = np.zeros(N)
    modes = []          # relevant modes
    sl = []
    
    # Two choices for limiting number of calculated modes
    # 1. Take all the modes for which S > S_limit
    # 2. Take 10 modes with the largest S
    # Option 1 ensures that all relevant modes are calculated, 
    # but option 2 is more consistent in calculation time. 
    # There might be large variance in number of relevant modes between different electronic states
    if 1: # 1.
        for i, f in enumerate(freq):
            S = get_huang_rhys(d[i], mu[i], f)
            if S > S_limit:
                sl.append((S, i))
    else: # 2.
        for i, f in enumerate(freq):
            S = get_huang_rhys(d[i], mu[i], f)
            sl.append((S, i))
        def comp(x): return x[0]
        sl = nlargest(10, sl, key=comp)
    
    ss = [x[0] for x in sl]
    modes = [x[1] for x in sl]
    print("\nMinimum required S: %14.10f" % min(ss))
    print("Number of modes in calculation: %i out of %i" % (len(modes), N))
    print("%10s %10s" % ("i", "S"))
    for S, i in sl:
        print("%10.i %10.8f" % (i, S))
    dat_I = []
    dat_E = []

    # itr_m is the list of combinations for ground state and itr_n for exited state
    # first argument in the combinations function is the most important number considering 
    # the script execution time. Anything above m=3 and n=4 not recommended since python isn't that fast.
    itr_m = combinations(m_lim, len(modes)) 
    itr_n = combinations(n_lim, len(modes))
   
    for mode in modes:
        if freq[mode] < 2:
            print("Negative mode as relevant")
    
    # Assuming that the calculation is from neutral ground state 
    # to ion and it's exited states (up in energy), then use the upper block.
    # If you want the transition down in energy, then uncomment the other block and comment out the first.
    # Blocks are identical in all other ways except in the way energy is calculated.
    #"""
    ##############################
    # Neutral -> Ion
    f = open(filename, "a")
    # Calculate intensity for all combinations of m and n
    for comb_m in itr_m:
        for i, mode in enumerate(modes):
            m[mode] = comb_m[i]
        # Get the energy of initial state
        e_m = get_vib_energy(freq, m, E0[0])
        for comb_n in itr_n:
            for i, mode in enumerate(modes):
                n[mode] = comb_n[i]
            I = intensity(d, freq, mu, m, n, T, S_limit, modes)
            # Filter out irrelevant transitions
            if I < 1e-3: continue
            # Get the energy of final state
            e_n = get_vib_energy(freq, n, E0[1])
            f.write("%10.6f  %10.6f  %s\n" % (I, e_n-e_m, c))
    #############################
    """
    #############################
    # Ion -> Neutral
    f = open(filename, "a")
    for comb_m in itr_m:
        for i, mode in enumerate(modes):
            m[mode] = comb_m[i]
        e_m = get_vib_energy(freq[0], m, E0[1])
        for comb_n in itr_n:
            for i, mode in enumerate(modes):
                n[mode] = comb_n[i]
            I = intensity(d, freq, mu, m, n, T, S_limit, modes)
            if I < 1e-3: continue
            e_n = get_vib_energy(freq[0], n, E0[0])
            f.write("%10.6f  %10.6f  %s\n" % (I, e_n-e_m, c))
    #############################
    #"""
    f.close()

def main(filename, tot_energy_file=""):
    
    rt = "vibrations/"
    folders_v = [rt+f+sub_fol for f in folders]
    print("Reading vibrational output")
    freq = get_frequencies(folders_v)                     
    mu = get_reduced_masses(folders_v)                    
    norm_modes = get_normal_modes(folders_v)              
    T = get_transformation_matrix(norm_modes[0])        
    pos = get_positions(folders_v)
    print("Done")

    """
    pos = []
    rt = "relaxations/"
    folders = ["neutral_0/", "ion_0/", "ion_1/", "ion_2/", "ion_3/", "ion_4/"]
    for fol in folders:
        f = open(rt+fol+"geometry.in.next_step")
        d = []
        for l in f:
            p = l.split()
            if p[0] == "atom":
                d.extend([float(t) for t in p[1:4]])
        pos.append(d)
        f.close()
    pos = [np.array(p) for p in pos]
    """

    # displacements
    disp = []
    for i in range(1, len(folders)):
        disp.append(pos[i]-pos[0])
    
    norm_disp = [np.dot(T, d) for d in disp]

    # displacement info
    for i, d in enumerate(disp):
        print("State %i displacement" % i)
        print("%14s %14.4f\n" % ("Magnitude", np.linalg.norm(d)))
        print("%14s %14s %14s" % ("Component", "Cartesian", "N. mode coord."))
        for j, v in enumerate(d):
            print("%14.i %14.4f %14.4f" % (j, v, norm_disp[i][j]))
        print("")

    # total energies
    E0 = []
    if tot_energy_file.split():
        f = open(tot_energy_file, "r")
        for l in f:
            E0.append(float(l))
        f.close()
    else:
        rt = "relaxations/"
        for fol in folders:
            f = open(rt+fol+"aims.out")
            for l in f:
                if " | Total energy of the DFT / Hartree-Fock s.c.f. calculation      : " in l:
                    E0.append(float(l.split()[11]))
                    break
            f.close()

    # clearing file
    f = open(filename, "w")
    f.write("%10s  %10s\n" % ("I","E (eV)"))
    f.close()
    
    # TODO: change naming to something sensible and define colors in plotting script
    c = ["blue", "red", "green", "indigo", "magenta", "aqua", "lime", "teal"]
   

    # writing to file
    for i, d in enumerate(disp):
        write_intensities(np.dot(T, d), freq[i+1], mu[i+1], [E0[0], E0[i+1]], 290, S_lim, c[i], filename)

    
    """
    # Using only the frequencies and normal modes of ground state
    # writing to file
    for i, d in enumerate(disp):
        write_intensities(np.dot(T, d), freq[0], mu[0], [E0[0], E0[i+1]], 290, S_lim, c[i], filename)
    """


if __name__ == "__main__":
    main(output_file, energy_file)



