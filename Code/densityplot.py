from matplotlib.pylab import *
from numpy import *

def n(r_val):
	return n0/(1 + (r_val/r0)**4)

n0 = 4.0
r0 = 10
r_val = linspace(0, 20, 1000)

f = open("density_profile.txt")

r = []
density = []

for line in f:
	words = line.split()
	r.append(float(words[0]))
	density.append(float(words[1]))

plot(r, density)
hold('on')
plot(r_val, n(r_val))
xlabel("radius R")
ylabel("Density")
show()