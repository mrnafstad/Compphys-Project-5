#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <time.h>
#include <vector>
#include "planet.h"
#include "solver.h"
#include "lib.h"
using namespace std;

// random numbers with gaussian distribution 
double gaussian_deviate(long * idum) 
{ 
	static int iset = 0; 
	static double gset; 
	double fac, rsq, v1, v2;

	if ( idum < 0) iset =0; 
	if (iset == 0) { 
		do { 
			v1 = 2.*ran2(idum) -1.0; 
			v2 = 2.*ran2(idum) -1.0; 
			rsq = v1*v1+v2*v2; 
		} while (rsq >= 1.0 || rsq == 0.); 
		fac = sqrt(-2.*log(rsq)/rsq);
		gset = v1*fac; 
		iset = 1; 
		return v2*fac; 
	} 	else { 
		iset =0; 
		return gset; 
	} 
} // end function for gaussian deviates


// A function that splits the ball of radius R0 into shells with shellsizes shell_size,
// and counts how many stars are found in each shell, writing it to file for plotting
void density_profile( solver allplanets, double shell_size, int N, double R0 )
{
	int thisshell, total = 0, bound = 0;
	double massthis;
	double r = 0, r_average1 = 0, r_average2 = 0, stddev1 = 0, stddev2 = 0;

	FILE *dens;
	char DensOut[30];
	snprintf(DensOut, sizeof(DensOut), "density_profile_%i.txt", N);

	dens = fopen(DensOut, "w+");

	//Loops until R0 is reached
	while ( r < R0) {
		thisshell = 0;
		massthis = 0;

		//loops over all planets 
		for ( int i = 0; i <= allplanets.total_planets; i++ ) {
			planet thisplanet = allplanets.all_planets[i];

			//If the planet isn't bound to the system, we simply skip it
			if(!thisplanet.bound) continue;

			double radius = 0, rad;

			//Calculates the radial coordinate of the planet
			for (int j = 0; j < 3; j++) radius += thisplanet.position[j]*thisplanet.position[j];
			rad  = sqrt(radius);

			//Adds together radi (only one time, hence e.g. r = 0.0)
			if (r == 0.0){
				r_average1 += rad;
				bound += 1;
				if(rad < R0 + shell_size) total += 1;
			}

			//Counting how many stars are found in the current shell
			if ( (rad > r) && (rad < r + shell_size) ) {
				thisshell += 1;
				massthis += thisplanet.mass;
				r_average2 += rad;
			}
		}

		//Calculates contributions to the standard deviations
		stddev1 += thisshell*pow(((r + shell_size/2) - r_average1/bound), 2);
		stddev2 += thisshell*pow(((r + shell_size/2) - r_average2/total), 2);
		//Writes the counted stars to file
		fprintf(dens, "%f %i \n", r, thisshell);
		r += shell_size;
	
	}

	//Prints out different calculations
	printf("\nStars within R0: %i \n", total);
	printf("Bound stars: %i \n", bound);
	printf("Average radius (stars inside R0): %.3f ly\n", r_average2/total);
	printf("Average radius (all bound stars): %.3f ly\n", r_average1/bound);
	printf("Standard deviation (starts inside R0): %.3f ly\n", sqrt(stddev2/total));
	printf("Standard deviation (all bound stars): %.3f ly\n", sqrt(stddev1/bound));
	fclose(dens);
}

// A function that checks if the virial theorem is fulfilled. It's only taking the bound planets into account
void virial ( solver allplanets )
{
	double kav = 0, uav = 0;
	int N = allplanets.total_planets, lost = 0;
	//Loops over all planets
	for ( int i = 0; i < N; i++ ) {
		planet thisplanet = allplanets.all_planets[i];
		//Checks whether the current planet is bound or not
		if ( thisplanet.bound ) {
			kav += thisplanet.kinetic;
			uav += thisplanet.potential;
		}
		else lost += 1;
	}
	double boundplanets = N - lost;

	printf("%f %f \n", (double)kav/boundplanets, -(double)uav/(double)boundplanets/2.0/2.0);
	//We divide by 2.0 twice, once for the virial theorem, and once since all potential energies is counted twice
}

int main ( int argc, char * argv[] )
{

	//NB N is now the number of initial bodies, dt is the number of integration points!
	int dim = 3, N;
	long idum = -1;

	double final_time, dt, R0, mass, R, x, y, z, theta, phi, vx = 0, vy = 0, vz = 0, shell_size;

	if(argc < 5){
		printf("Bad usage, cml arguments: N R0 final_time \n");
		return 0;
	}
	else{
		N = atoi(argv[1]);
		R0 = atof(argv[2]);
		final_time = atof(argv[3]);
		dt = atof(argv[4]);
	}
	//For use in virial()
	shell_size = R0/(double)100;

	/*---------------------------
	IMPORTANT!
	The following bool statements are deciding how the program runs:
	-bool energy: if true it calculates and prints out the total energies a certain time intervals to terminal
	-bool smoothing: if true it uses a modified newtonian force found in solver.cpp
	-bool constantMass: if true a constant mass M0 defined below is set to the total mass of the system, regardless of N
						The average mass per particle is then calculated from this value
	*/
	bool energy = true;
	bool smoothing = true;
	bool constantMass = true;

	// Only used if constantMass = true
	double M0 = 1000;  // solar masses 
	double mu = M0/N;  // average mass per particle

	// creates object of the solver.cpp file
	solver system_VV(R0);

	// Loops over all particles
	for ( int i = 0; i < N; i++) {
		// Calculates (uniformly) random postions of the particles within a sphere of radius R0
		R = ran1(&idum)*R0;
		theta = ran1(&idum)*2*M_PI;
		phi = ran1(&idum)*M_PI;
		x = R*cos(theta)*sin(phi);
		y = R*sin(theta)*sin(phi);
		z = R*cos(phi);


		if(constantMass){
			// Random mass around mu solar masses
			// 100/N adjusting stddev so that mass is never negative
			mass = 100/(double)N*gaussian_deviate(&idum) + mu;
		}
		else{
			// Random mass around 10 solar masses with stddev 1
			mass = gaussian_deviate(&idum) + 10;
		}

		// Adding the planets to the system
		planet thistest( mass, x, y, z, vx, vy, vz );
		system_VV.add( thistest );
	}

	printf("Total stars: %i\n", system_VV.total_planets);

	system_VV.GravitationalConstant(); // Setting the new gravitational constant found in solver.cpp
	system_VV.velVerlet( dim, dt, final_time, N, energy, smoothing); // Calculates the motion of the system using the Velocity Verlet method
	density_profile( system_VV, shell_size, N, R0); // Finds the density profile of the system at the last time step
	virial(system_VV); //Checks if the virial theorem is fulfilled at the last timestep

	return 0;
}
