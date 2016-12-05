from matplotlib.pylab import *
from numpy import *

f = open("VerletTest.txt")

lines = f.readlines()

words = lines[0].split()
num =  len(words) # Numbers per line
stars = (num -1)/3       # Three coordinates per galaxies


first_line = lines[:1][0].split()
last_line = lines[-1].split()

last_line = array(last_line[1:])

x_first = []; x_last = []
y_first = []; y_last = []

for i in range(stars):
	x_first.append(float(first_line[3*i +1]))
	y_first.append(float(first_line[3*i +2]))
	x_last.append(float(last_line[3*i +1]))
	y_last.append(float(last_line[3*i +2]))

x_first = array(x_first); x_last = array(x_last)
y_first = array(y_first); y_last = array(y_last)

plot(x_first, y_first, 'o')
show()
plot(x_last, y_last, 'o')
xlim([-50, 50])
ylim([-50, 50])
show()

f.close()