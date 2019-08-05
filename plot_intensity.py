import matplotlib.pyplot as plt
import numpy as np

# xlim, ylim
mi = 9.0
ma = 15.0
# sigma for gaussian
sig = 0.05
# line resolution
dx = 0.001
x = np.arange(mi, ma, dx)

# sigmas for plot_gauss()
sigl = [0.05, 0.08, 0.15, 0.3]

# colors for states when plotting
cl = ["blue","red","green","indigo","magenta","aqua","lime","teal"]

# gaussian function
def gaussian(x, mu, sig):
    return 1/(sig*np.sqrt(2*np.pi)) * np.exp(-np.power(x - mu, 2) / (2 * np.power(sig, 2)))

# some figures to begin with
fig1 = plt.figure(figsize=(13, 11))
#fig2 = plt.figure(figsize=(13, 11))
fig3 = plt.figure(figsize=(13, 11))

# list for plot_zz()
tot_gauss = []

# Plots the delta peaks by themselves and with gaussian broadening
# parameters: data file, figure, title
def plot_data(filename, fig, title):
    peaks_only = 1

    ax = fig.add_subplot(211)
    ax12 = fig.add_subplot(212, sharex=ax)
    fig.subplots_adjust(hspace=0.1)
    dat_I = []
    dat_E = []
    col = []
    print("Reading file")
    f = open(filename, "r")
    f.next()

    zz_peaks = []   # (I, E)
    dic = {}
    ind = 0

    for l in f:
        p = l.split()
        dat_I.append(float(p[0]))
        dat_E.append(abs(float(p[1])))
        c = int(p[2])
        if c not in col:
            zz_peaks.append((float(p[0]), abs(float(p[1]))))
            dic[c] = ind
            ind += 1
        col.append(c)
    f.close()
    
    print("Reading done")
    n = len(dic)

    inv_dic = {v: k for k, v in dic.iteritems()}

    y = [np.zeros(len(x)) for _ in range(n)]
    y0 = [np.zeros(len(x)) for _ in range(n)]

    print("Starting gaussian broadening")
    for j, (i, e, c) in enumerate(zip(dat_I, dat_E, col)):
        if j % 10 == 0:
            print(j)
        for j, v in enumerate(x):
            y[dic[c]][j] += i*gaussian(v, e, sig)
    """
    for i in range(n):
        for j, v in enumerate(x):
            y0[i][j] += zz_peaks[i][0]*gaussian(v, zz_peaks[i][1], sig)
    """
    print("Gaussian broadening ready")

    max_all = [np.argmax(y[i]) for i in range(n)]
    max_0 = [np.argmax(y0[i]) for i in range(n)]

    print("Plotting")
    """
    for i in range(n):
        ax.plot([mi+max_0[i]*dx, mi+max_0[i]*dx], [0, y0[i][max_0[i]]], "--", color=inv_dic[i])
        ax.plot([mi+max_all[i]*dx, mi+max_all[i]*dx], [0, y[i][max_all[i]]], "--", color=inv_dic[i], label="Peak diff. %i: %4.3f eV" % (i, abs(x[max_0[i]]-x[max_all[i]])))
    """
    for i in range(len(dat_I)):
        if peaks_only:
            ax12.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)
        ax.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)

    for i, c in enumerate(set(col)):
        ax.plot(x, y[dic[c]], "-", color=cl[c], label="")
        #ax.plot(x, y0[dic[c]], "-.", lw=1.5, color=cl[c])

    y_all = np.zeros(len(x))
    for i, v in enumerate(x):
        y_all[i] = sum([y[j][i] for j in range(n)])

    ax.plot(x, y_all, "-", color="black")
    
    tot_gauss.append(y_all)

    ax.set_xlim(mi, ma)
    #ax.set_ylim(0, 3)
    ax.set_title(title)
    #ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Relative intensity")
    
    ax12.set_xlim(mi, ma)
    #ax12.set_ylim(0, 0.5)
    ax12.set_xlabel("Energy (eV)")
    ax12.set_ylabel("Relative intensity")
    
    ax.legend()

# Plots several spectrum on the same figure. 
# Not maintained.
def plot_zz(fig):
    pbe = []
    dci = []
    b3l = []

    f = open("DCI.dat", "r")
    for l in f:
        dci.append(float(l))
    f.close()
    
    f = open("pbe_zz.dat", "r")
    for l in f:
        pbe.append(float(l))
    f.close()

    f = open("b3lyp.dat", "r")
    for l in f:
        b3l.append(float(l))
    f.close()
    
    y_pbe = np.zeros(len(x))
    for e in pbe:
        for i, v in enumerate(x):
            y_pbe[i] += gaussian(v, e, sig)
    tot_gauss.append(y_pbe)

    y_dci = np.zeros(len(x))
    for e in dci:
        for i, v in enumerate(x):
            y_dci[i] += gaussian(v, e, sig)
    tot_gauss.append(y_dci)
    
    y_b3l = np.zeros(len(x))
    for e in b3l:
        for i, v in enumerate(x):
            y_b3l[i] += gaussian(v, e, sig)
    tot_gauss.append(y_b3l)
    
    fig3.subplots_adjust(hspace=0)

    ax = fig3.add_subplot(111)
    ax.set_xlim(mi, ma)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.set_yticks([])

    ax_0 = fig.add_subplot(611)
    ax_0.xaxis.set_ticklabels([])
    ax_0.set_yticks([])
    ax_0.set_xlim(mi, ma)
    
    ax_1 = fig.add_subplot(612, sharex=ax_0)
    ax_1.xaxis.set_ticklabels([])
    ax_1.set_yticks([])
    
    ax_2 = fig.add_subplot(615, sharex=ax_0)
    ax_2.xaxis.set_ticklabels([])
    ax_2.set_yticks([])

    ax_3 = fig.add_subplot(616, sharex=ax_0)
    ax_3.set_yticks([])

    ax_4 = fig.add_subplot(613, sharex=ax_0)
    ax_4.xaxis.set_ticklabels([])
    ax_4.set_yticks([])
    
    ax_5 = fig.add_subplot(614, sharex=ax_0)
    ax_5.xaxis.set_ticklabels([])
    ax_5.set_yticks([])

    ax_0.plot(x, tot_gauss[0], "-", c="b", label="PBE, vib")
    ax_0.set_ylim(0, max(tot_gauss[0])*1.1)
    
    ax_1.plot(x, tot_gauss[3], "-", c="b", label="PBE, 0-0")
    ax_1.set_ylim(0, max(tot_gauss[3])*1.1)
    
    ax_4.plot(x, tot_gauss[1], "-", c="b", label="B3LYP, vib")
    ax_4.set_ylim(0, max(tot_gauss[1])*1.1)
    
    ax_5.plot(x, tot_gauss[5], "-", c="b", label="B3LYP, 0-0")
    ax_5.set_ylim(0, max(tot_gauss[5])*1.1)
    
    ax_2.plot(x, tot_gauss[2], "-", c="b", label="DCI & PBE, vib")
    ax_2.set_ylim(0, max(tot_gauss[2])*1.1)
    
    ax_3.plot(x, tot_gauss[4], "-", c="b", label="DCI, 0-0")
    ax_3.set_ylim(0, max(tot_gauss[4])*1.1)

    ax_0.legend()
    ax_1.legend()
    ax_2.legend()
    ax_3.legend()
    ax_4.legend()
    ax_5.legend()

    ax.set_title("Neutral -> Ion")
    ax.set_ylabel("Intensity")
    ax.set_xlabel("Energy (eV)")
    ax.legend()

# Same as plotting() but places of 0-0 peaks are shifted accordin to shift_file.
def plot_data_shift(filename, shift_file, fig, title):
    peaks_only = 1

    ax = fig.add_subplot(211)
    ax12 = fig.add_subplot(212, sharex=ax)
    fig.subplots_adjust(hspace=0.1)
    dat_I = []
    dat_E = []
    col = []
    print("Reading file")
    f = open(filename, "r")
    f.next()

    zz_peaks = []   # (I, E)
    dic = {}
    ind = 0

    for l in f:
        p = l.split()
        dat_I.append(float(p[0]))
        dat_E.append(abs(float(p[1])))
        c = int(p[2])
        if c not in col:
            zz_peaks.append((float(p[0]), abs(float(p[1]))))
            dic[c] = ind
            ind += 1
        col.append(c)
    f.close()

    dci = []

    f = open(shift_file, "r")
    for l in f:
        dci.append(float(l))
    f.close()

    for i in range(len(dat_E)):
        dat_E[i] -= zz_peaks[dic[col[i]]][1]
        dat_E[i] += dci[dic[col[i]]+1]

    print("Reading done")
    n = len(dic)

    inv_dic = {v: k for k, v in dic.iteritems()}

    y = [np.zeros(len(x)) for _ in range(n)]
    y0 = [np.zeros(len(x)) for _ in range(n)]

    print("Starting gaussian broadening")
    for j, (i, e, c) in enumerate(zip(dat_I, dat_E, col)):
        if j % 10 == 0:
            print(j)
        for j, v in enumerate(x):
            y[dic[c]][j] += i*gaussian(v, e, sig)
    """
    for i in range(n):
        for j, v in enumerate(x):
            y0[i][j] += zz_peaks[i][0]*gaussian(v, zz_peaks[i][1], sig)
    """
    print("Gaussian broadening ready")

    max_all = [np.argmax(y[i]) for i in range(n)]
    max_0 = [np.argmax(y0[i]) for i in range(n)]

    print("Plotting")
    """
    for i in range(n):
        ax.plot([mi+max_0[i]*dx, mi+max_0[i]*dx], [0, y0[i][max_0[i]]], "--", color=inv_dic[i])
        ax.plot([mi+max_all[i]*dx, mi+max_all[i]*dx], [0, y[i][max_all[i]]], "--", color=inv_dic[i], label="Peak diff. %i: %4.3f eV" % (i, abs(x[max_0[i]]-x[max_all[i]])))
    """
    for i in range(len(dat_I)):
        if peaks_only:
            ax12.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)
        ax.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)

    for i, c in enumerate(set(col)):
        ax.plot(x, y[dic[c]], "-", color=cl[c], label="")
        #ax.plot(x, y0[dic[c]], "-.", lw=1.5, color=cl[c])

    y_all = np.zeros(len(x))
    for i, v in enumerate(x):
        y_all[i] = sum([y[j][i] for j in range(n)])

    ax.plot(x, y_all, "-", color="black")
    
    tot_gauss.append(y_all)

    ax.set_xlim(mi, ma)
    #ax.set_ylim(0, 3)
    ax.set_title(title)
    #ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Relative intensity")
    
    ax12.set_xlim(mi, ma)
    #ax12.set_ylim(0, 0.5)
    ax12.set_xlabel("Energy (eV)")
    ax12.set_ylabel("Relative intensity")
    
    ax.legend()



def plot_shift_vib(filename, shift_file, ax, title):

    dat_I = []
    dat_E = []
    col = []
    print("Reading file")
    f = open(filename, "r")
    f.next()

    zz_peaks = []
    dic = {}
    ind = 0

    for l in f:
        p = l.split()
        dat_I.append(float(p[0]))
        dat_E.append(abs(float(p[1])))
        c = p[2]
        if c not in col:
            zz_peaks.append((float(p[0]), abs(float(p[1]))))
            dic[c] = ind
            ind += 1
        col.append(c)
    f.close()

    dci = []

    f = open(shift_file, "r")
    for l in f:
        dci.append(float(l))
    f.close()

    for i in range(len(dat_E)):
        dat_E[i] -= zz_peaks[dic[col[i]]][1]
        dat_E[i] += dci[dic[col[i]]+1]
    
    print("Reading done")
    n = len(dic)

    inv_dic = {v: k for k, v in dic.iteritems()}

    y = [np.zeros(len(x)) for _ in range(n)]

    print("Starting gaussian broadening")
    for j, (i, e, c) in enumerate(zip(dat_I, dat_E, col)):
        if j % 10 == 0:
            print(j)
        for j, v in enumerate(x):
            y[dic[c]][j] += i*gaussian(v, e, sig)

    print("Gaussian broadening ready")

    max_all = [np.argmax(y[i]) for i in range(n)]

    print("Plotting")

    for i in range(len(dat_I)):
        ax22.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)
        ax.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], color=cl[col[i]], lw=2.0)

    for c in set(col):
        ax.plot(x, y[dic[c]], "-", color=c)

    y_all = np.zeros(len(x))
    for i, v in enumerate(x):
        y_all[i] = sum([y[j][i] for j in range(n)])
    
    tot_gauss.append(y_all)
    ax.plot(x, y_all, "-", color="black")

    ax.set_title(title)

    ax.set_xlim(mi, ma)
    #ax.set_ylim(0, 3)
    #ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Relative intensity")
    
    ax22.set_xlim(mi, ma)
    #ax22.set_ylim(0, 0.5)
    ax22.set_xlabel("Energy (eV)")
    ax22.set_ylabel("Relative intensity")
    
    ax.legend()

# plots the state alone and calculates the shift from 0-0 for different sigma
# parameters: data file, state, subplot
def plot_gauss(filename, tr, ax):
    dat_I = []
    dat_E = []
    print("Reading file")
    f = open(filename, "r")
    f.next()
    
    col = ["r","g","c","m"]

    zz_peaks = []   # (I, E)

    for l in f:
        p = l.split()
        if int(p[2]) == tr:
            dat_I.append(float(p[0]))
            dat_E.append(abs(float(p[1])))
    f.close()
    zz_peaks = (dat_I[0], dat_E[0])
    
    print("Reading done")
    
    n = 4

    y = [np.zeros(len(x)) for _ in range(n)]
    y0 = np.zeros(len(x))

    print("Starting gaussian broadening")
    for k in range(n):
        for j, (i, e) in enumerate(zip(dat_I, dat_E)):
            if j % 10 == 0:
                print(j)
            for j, v in enumerate(x):
                y[k][j] += i*gaussian(v, e, sigl[k])
    
    for j, v in enumerate(x):
        y0[j] += zz_peaks[0]*gaussian(v, zz_peaks[1], sig)
    
    print("Gaussian broadening ready")

    max_all = [np.argmax(y[i]) for i in range(n)]
    max_0 = np.argmax(y0)

    print("Plotting")
    
    ax.plot([mi+max_0*dx, mi+max_0*dx], [0, y[0][max_all[0]]], "--")
    for i in range(n):
        ax.plot([mi+max_all[i]*dx, mi+max_all[i]*dx], [0, y[i][max_all[i]]], "--", c=col[i])
    for i in range(len(dat_I)):
        ax.plot([dat_E[i], dat_E[i]], [0, dat_I[i]], c="b", lw=2.0)

    lab = ["a","b","c","d"]

    for i in range(n):
        #ax.plot(x, y[i], "-", c=col[i], label=r'$d_%s = %4.3f \: \mathrm{eV}$' % (lab[i], abs(x[max_0]-x[max_all[i]])))
        ax.plot(x, y[i], "-", c=col[i], label='d_%s = %4.3f eV' % (lab[i], abs(x[max_0]-x[max_all[i]])))

    ax.set_xlim(mi+max_0*dx-0.6, mi+max_0*dx+1.4)
    #ax.set_ylim(0, 3)
    ax.set_yticks([])
    #ax.set_xticks([])
    ax.legend()

def plot_peak_shift(filename, states, fig):

    ax9 = fig.add_subplot(111)
    ax9.set_xlim(mi, ma)
    ax9.spines['top'].set_color('none')
    ax9.spines['bottom'].set_color('none')
    ax9.spines['left'].set_color('none')
    ax9.spines['right'].set_color('none')
    ax9.set_yticks([])
    ax9.set_xticks([])
    ax9.set_xlabel("Energy (eV)")

    #ax9.set_title("Peak shift")
    # Latex in pyplot might not work
    ax9.set_title(r'$\mathrm{Sigma \  comparison} \quad \sigma_a=0.05 \quad \sigma_b=0.08 \quad \sigma_c=0.15 \quad \sigma_d=0.3$')
    ax9.xaxis.set_label_coords(0.5, -0.04)
    fig.subplots_adjust(hspace=0.1)
    
    for s in states:
        sp = str(len(states))+"1"+str(s+1)
        plot_gauss(filename, s, fig.add_subplot(int(sp)))



if __name__ == "__main__":
   
    # Plots peaks from a file to figure
    plot_data("intensity.dat", fig1, "")
   
    # Plots states on their own and calculates the shift from 0-0 transition.
    plot_peak_shift("intensity.dat", [0, 1, 2], fig3)

   
    # Saving figure
    #fig1.savefig("results/neutral-ion_dci.png")
    
    
    plt.show()

