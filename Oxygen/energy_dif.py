import numpy as np

f = open("GW/aims.out", "r")

n = 9
eigval = []

for l in f:
    if "GW quasi-particle energy levels" in l:
        for _ in range(6): next(f)
        eigval = []
        for _ in range(n):
            eigval.append(float(next(f).split()[6]))
f.close()

print(eigval)

neu = 2*sum(eigval)
states = range(7)

f = open("zz_energies.dat", "w")
f.write("%10.4f\n" % 0)
for s in states:
    enr = 0
    for i, val in enumerate(eigval[::-1]):
        if i == s:
            enr += val
        else:
            enr += 2*val
    f.write("%10.4f\n" % (enr-neu))

f.close()
