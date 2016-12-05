from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D
from numpy import *

f = open("VerletTest.txt")

words = f.readline().split()
num =  len(words) # Numbers per line
galaxies = (num -1)/3       # Three coordinates per galaxies

x = zeros(galaxies); y = zeros(galaxies) ; z = zeros(galaxies)

fig = figure()
ax = fig.add_subplot(111, projection='3d')

ion()
show()

for line in f:
	words = line.split()
	t = float(words[0])
	for i in range(galaxies):
		x[i] = float(words[3*i +1])
		y[i] = float(words[3*i +2])
		z[i] = float(words[3*i +2])

	ax.cla()
	ax.plot(x,y, z, 'o')
	ax.set_xlim([-20, 20])
	ax.set_ylim([-20, 20])
	ax.set_zlim([-20, 20])
	ax.legend(["t = %.2f" % t])
	pause(0.01)

f.close()
