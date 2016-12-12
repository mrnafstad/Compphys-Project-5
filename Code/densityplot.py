from matplotlib.pylab import *
from numpy import *

def n(r_val, n0, r0):
	return n0/(1 + (r_val/r0)**4)

n0_1 = 3;  n0_2 = 15;  n0_3 = 25
r0_1 = 10;   r0_2 = 7 ;  r0_3 = 6.5

r_val = linspace(0, 20, 1000)

f1 = open("density_profile_100.txt")
f2 = open("density_profile_500.txt")
f3 = open("density_profile_1000.txt")

r1 = []; r2 = []; r3 = []
density1 = []; density2 = []; density3 = []

for line in f1:
	words = line.split()
	r1.append(float(words[0]))
	density1.append(float(words[1]))

for line in f2:
	words = line.split()
	r2.append(float(words[0]))
	density2.append(float(words[1]))

for line in f3:
	words = line.split()
	r3.append(float(words[0]))
	density3.append(float(words[1]))

r1 = array(r1); density1 = array(density1)
r2 = array(r2); density2 = array(density2)
r3 = array(r3); density3 = array(density3)

plot(r1, density1/len(r1))
hold('on')
plot(r2, density2/len(r2))
hold('on')
plot(r3, density3/len(r3))
#plot(r_val, n(r_val)/len(r))
xlabel("R")
ylabel("Density")
legend(["N = 100", "N = 500", "N = 1000"])
show()

plot(r1, density1)#/len(r1))
hold('on')
plot(r_val, n(r_val, n0_1, r0_1))
xlabel("R")
ylabel("Density")
show()

plot(r2, density2)#/len(r2))
hold('on')
plot(r_val, n(r_val, n0_2, r0_2))
xlabel("R")
ylabel("Density")
show()

plot(r3, density3)#/len(r3))
hold('on')
plot(r_val, n(r_val, n0_3, r0_3))
xlabel("R")
ylabel("Density")
show()