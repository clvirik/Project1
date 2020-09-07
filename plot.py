import matplotlib.pyplot as plt
import glob
import numpy as np

txtfiles = []
for file in glob.glob("N*.txt"):
    txtfiles.append(file)

plt.figure(1)

for file in txtfiles:
    v = []; u = []; x = []

    infile = open(file, 'r')
    for line in infile:
        words = line.split()
        x.append(float(words[0]))
        v.append(float(words[1]))
        u.append(float(words[2]))
    infile.close()

    plt.plot(x, v, label='v(x)')

plt.plot(x, u, label='u(x)')
plt.xlabel('x'); plt.ylabel('y')
plt.legend()
plt.savefig("plot.png")

infile = open("error.txt", 'r')
N = []; error = []
for line in infile:
    words = line.split()
    N.append(float(words[0]))
    error.append(float(words[1]))

plt.figure(2)

plt.plot(np.log10(np.asarray(N)), error, label=r'$\epsilon (N)$')
plt.xlabel(r'$\log_{10}(N)$'); plt.ylabel(r'$\log_{10}(\epsilon)$')
plt.legend()
plt.savefig("error.png")