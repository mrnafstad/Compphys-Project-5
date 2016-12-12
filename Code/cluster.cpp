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

void density_profile( solver allplanets, double shell_size, int N, double R0 )
{
	int thisshell, total = 0, bound = 0;
	double massthis;
	double r = 0, r_average1 = 0, r_average2 = 0, stddev1 = 0, stddev2 = 0;

	FILE *dens;
	char DensOut[30];
	snprintf(DensOut, sizeof(DensOut), "density_profile_%i.txt", N);

	dens = fopen(DensOut, "w+");

	while ( r < R0) {
		thisshell = 0;
		massthis = 0;
		//std::cout << allplanets.total_planets << std::endl;

		for ( int i = 0; i <= allplanets.total_planets; i++ ) {
			planet thisplanet = allplanets.all_planets[i];

			if(!thisplanet.bound) continue;

			double radius = 0, rad;

			for (int j = 0; j < 3; j++) radius += thisplanet.position[j]*thisplanet.position[j];
			rad  = sqrt(radius);

			//Adds together radi (only one time, hence e.g. r = 0.0)
			if (r == 0.0){
				r_average1 += rad;
				bound += 1;
				if(rad < R0 + shell_size) total += 1;
			}

			if ( (rad > r) && (rad < r + shell_size) ) {
				thisshell += 1;
				massthis += thisplanet.mass;
				r_average2 += rad;
			}
		}

		stddev1 += thisshell*pow(((r + shell_size/2) - r_average1/bound), 2);
		stddev2 += thisshell*pow(((r + shell_size/2) - r_average2/total), 2);
		fprintf(dens, "%f %i \n", r, thisshell);
		r += shell_size;
	
	}

	printf("\nStars within R0: %i \n", total);
	printf("Bound stars: %i \n", bound);
	printf("Average radius (stars inside R0): %.3f ly\n", r_average2/total);
	printf("Average radius (all bound stars): %.3f ly\n", r_average1/bound);
	printf("Standard deviation (starts inside R0): %.3f ly\n", sqrt(stddev2/total));
	printf("Standard deviation (all bound stars): %.3f ly\n", sqrt(stddev1/bound));
	fclose(dens);
}

void virial ( solver allplanets )
{
	double kav = 0, uav = 0;
	int N = allplanets.total_planets, lost = 0;
	for ( int i = 0; i < N; i++ ) {
		planet thisplanet = allplanets.all_planets[i];
		if ( thisplanet.bound ) {
			kav += thisplanet.kinetic;
			uav += thisplanet.potential;
		}
		else lost += 1;
	}
	double boundplanets = N - lost;
	std::cout << lost << std::endl;

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
	shell_size = R0/(double)100;

	bool energy = true;
<<<<<<< HEAD
	bool smoothing = true;
	bool constantMass = true;

	double M0 = 1000;  // solar masses 
	double mu = M0/N;  // average mass per particle
=======
	bool smoothing = false;
>>>>>>> eb4e2444095254ab47471ac2a9f0ae4299334cc4

	solver system_VV(R0);

	for ( int i = 0; i < N; i++) {
		R = ran1(&idum)*R0;
		theta = ran1(&idum)*2*M_PI;
		phi = ran1(&idum)*M_PI;
		x = R*cos(theta)*sin(phi);
		y = R*sin(theta)*sin(phi);
		z = R*cos(phi);

		if(constantMass){
			// 100/N adjusting stddev so that mass is never negative
			mass = 100/(double)N*gaussian_deviate(&idum) + mu;
		}
		else{
			mass = gaussian_deviate(&idum) + 10;
		}

		planet thistest( mass, x, y, z, vx, vy, vz );
		system_VV.add( thistest );
	}

	printf("Total stars: %i\n", system_VV.total_planets);


	system_VV.GravitationalConstant();
	system_VV.velVerlet( dim, dt, final_time, N, energy, smoothing);
	//printf("%lf\n", system_VV.all_planets[0].position[0]);
	density_profile( system_VV, shell_size, N, R0);
	virial(system_VV);
	return 0;
}