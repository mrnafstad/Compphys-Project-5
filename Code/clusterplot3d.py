from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D
from numpy import *


f = open("VerletTest.txt")
#f = open("clusterwithsmooth.txt")
#f = open("clusterwithout.txt")

words = f.readline().split()
num =  len(words) # Numbers per line
galaxies = (num -1)/3       # Three coordinates per galaxies

x = zeros(galaxies); y = zeros(galaxies) ; z = zeros(galaxies)


fig = figure()
ax = fig.add_subplot(111, projection='3d')

ion()
show()

# Plots the positions of the particles (3D) continously, creating an animation
for line in f:
	words = line.split()
	t = float(words[0])
	for i in range(galaxies):
		x[i] = float(words[3*i +1])
		y[i] = float(words[3*i +2])
		z[i] = float(words[3*i +3])

	ax.cla()
	ax.plot(x,y, z, 'o')
	ax.set_xlim([-20, 20])
	ax.set_ylim([-20, 20])
	ax.set_zlim([-20, 20])
	ax.set_xlabel("x (ly)")
	ax.set_ylabel("y (ly)")
	ax.set_zlabel("z (ly)")
	ax.legend(["t = %.2f" % t])
	if t < 0.02:
		pause(10)
	elif t == 1:
		pause(10)
	elif t == 3:
		pause(10)
	else:
		pause(0.001)




f.close()
