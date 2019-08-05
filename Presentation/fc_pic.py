import matplotlib.pyplot as plt
import numpy as np
import scipy.special as spec
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection


def func(R, n, l, a):
    N = np.sqrt(np.sqrt(a/3) * 1/(2**n * spec.factorial(n)))
    return N*0.5*hermite_polynomial(n, np.sqrt(a)*(R-l))*np.exp(-a*(R-l)**2 / 2)

def hermite_polynomial(n, x):
    p = spec.hermite(n)
    return p(x)

def plot_tr(n):
    
    k = 50

    fig = plt.figure(figsize=(9, 12))
    xx = np.arange(0.2, 1.8, 0.001)

    ax1 = fig.add_subplot(111)
    p = np.poly1d([k, 0, 0])
    ax1.plot(xx, p(xx-1), "-", c="teal")

    xx1 = np.array([x+0.3 for x in xx])
    p = np.poly1d([k, 0, 40])
    ax1.plot(xx1, p(xx1-1.3), "-", c="teal")

    for i in range(9):
        xf = np.arange(0.8-0.2*np.sqrt(i), 1.2+0.2*np.sqrt(i), 0.001)
        if i == 0:
            ax1.plot(xf, 1+2.7*i+(func(xf, i, 1.0, 40)), c="r")
            continue
        ax1.plot(xf, 1+2.7*i+(func(xf, i, 1.0, 40)), c="g")

    for i in range(9):
        xf = np.arange(1.1-0.2*np.sqrt(i), 1.5+0.2*np.sqrt(i), 0.001)
        if i == n:
            ax1.plot(xf, 41+2.7*i+(func(xf, i, 1.3, 40)), c="r")
            continue
        ax1.plot(xf, 41+2.7*i+(func(xf, i, 1.3, 40)), c="g")
   
    ax1.arrow(1, 1, 0, 40+2.7*n, color="indigo", length_includes_head=True,  head_length=2, head_width=0.04)

    #ax1.plot([0, 0],[-10, 0],"--", c="k")
    #ax1.plot([1.0, 1.0],[-10, 0],"--", c="k")
    #ax1.plot([1.3, 1.3],[-10, 0],"--", c="k")

    ax1.set_xlim(-0.5, 2.5)
    ax1.set_ylim(-10, 90)

    ax1.set_ylabel("Energy", fontsize=24)
    ax1.set_xlabel("Displacement", fontsize=24)
    
    ax1.set_xticks([])
    ax1.set_yticks([])

    fig.tight_layout()
    
    
    fig.savefig("fc/tr_%i.png" % n)

def plot_atoms():

    k = 50

    pat = [Ellipse((1, -5), 0.18, 4), Ellipse((1.3, -5), 0.18, 4), Ellipse((0, -5), 0.18, 4)]

    fig = plt.figure(figsize=(9, 12))
    xx = np.arange(0.2, 1.8, 0.001)

    ax1 = fig.add_subplot(111)
    p = np.poly1d([k, 0, 0])
    ax1.plot(xx, p(xx-1), "-", c="teal")

    xx1 = np.array([x+0.3 for x in xx])
    p = np.poly1d([k, 0, 40])
    ax1.plot(xx1, p(xx1-1.3), "-", c="teal")
    
    for i in range(9):
        xf = np.arange(0.8-0.2*np.sqrt(i), 1.2+0.2*np.sqrt(i), 0.001)
        if i == 0:
            ax1.plot(xf, 1+2.7*i+(func(xf, i, 1.0, 40)), c="r")
            continue
        ax1.plot(xf, 1+2.7*i+(func(xf, i, 1.0, 40)), c="g")
    
    for i in range(9):
        xf = np.arange(1.1-0.2*np.sqrt(i), 1.5+0.2*np.sqrt(i), 0.001)
        if i == 2:
            ax1.plot(xf, 41+2.7*i+(func(xf, i, 1.3, 40)), c="r")
            continue
        ax1.plot(xf, 41+2.7*i+(func(xf, i, 1.3, 40)), c="g")
    
    ax1.arrow(1, 1, 0, 40+2.7*2, color="indigo", length_includes_head=True,  head_length=2, head_width=0.04)

    p = PatchCollection(pat[0::1], alpha=0.9, color="red")
    ax1.add_collection(p)

    ax1.set_xlim(-0.5, 2.5)
    ax1.set_ylim(-10, 90)

    ax1.set_ylabel("Energy", fontsize=24)
    ax1.set_xlabel("Distance", fontsize=24)
    
    fig.tight_layout()
    ax1.tick_params(axis='both', which='major', labelsize=22)
    ax1.tick_params(axis='both', which='minor', labelsize=22)

    ax1.set_xticks([])
    ax1.set_yticks([])

    fig.savefig("fc_build/step7.png")


def main():
#    plot_atoms()

    for i in range(9):
        plot_tr(i)

if __name__ == "__main__":
    main()

