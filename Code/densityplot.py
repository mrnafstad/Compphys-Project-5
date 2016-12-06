from matplotlib.pylab import *
from numpy import *

f = open("density_profile.txt")

r = []
density = []

for line in f:
	words = line.split()
	r.append(float(words[0]))
	density.append(float(words[1]))

plot(r, density)
xlabel("radius R")
ylabel("Density")
show()