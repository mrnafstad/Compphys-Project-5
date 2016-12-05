from matplotlib.pylab import *
from numpy import *

f = open("VerletTest.txt")

#lines = f.readlines()

words = f.readline().split()
num =  len(words) # Numbers per line
galaxies = (num -1)/3       # Three coordinates per galaxies

x = zeros(galaxies); y = zeros(galaxies)

ion()
show()

for line in f:
	words = line.split()
	t = float(words[0])
	for i in range(galaxies):
		x[i] = float(words[3*i +1])
		y[i] = float(words[3*i +2])
	clf()
	plot(x,y, 'o')
	legend(["t = %.2f" %t])
	xlim([-20, 20])
	ylim([-20, 20])
	pause(0.01)

f.close()