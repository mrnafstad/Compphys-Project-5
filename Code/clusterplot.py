from matplotlib.pylab import *
from numpy import *

f = open("VerletTest.txt")

words = f.readline().split()
num =  len(words) # Numbers per line
galaxies = (num -1)/3       # Three coordinates per galaxies

x = [] ; y = [];

for i in range(galaxies):
	x.append(words[3*i +1])
	y.append(words[3*i +2])

x = array(x) ; y = array(y);

plot(x, y, 'o')
show()