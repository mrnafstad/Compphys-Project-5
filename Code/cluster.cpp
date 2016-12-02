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



int main ( int argc, char * argv[] )
{
	bool stationary = false;
	bool energy = false;
	bool relativity = false;
	bool MercPeri = false;
	//NB N is now the number of initial bodies, dt is the number of integration points!
	int dim = 3, dt = atoi(argv[1]), N = atoi(argv[2]);
	long idum = -1;

	double final_time =  atof(argv[3]), R0 = atof(argv[4]), M0 = 1.0;
	double mass, R, x, y, z, theta, phi, vx = 0, vy = 0, vz = 0;

	vector<planet> cluster;

	solver system_VV;

	for ( int i = 0; i <= N; i++) {
		R = ran1(&idum)*R0;
		theta = ran1(&idum)*2*M_PI;
		phi = ran1(&idum)*M_PI;
		x = R*cos(theta)*sin(phi);
		y = R*sin(theta)*sin(phi);
		z = R*cos(phi);

		mass = ran1(&idum);

		planet thistest( mass, x, y, z, vx, vy, vz );
		cluster.push_back( thistest );
		system_VV.add( thistest );
	}

	system_VV.velVerlet( dim, dt, final_time, energy, stationary, relativity, MercPeri);

	
	return 0;
}