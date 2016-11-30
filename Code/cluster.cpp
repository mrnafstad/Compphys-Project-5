#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <time.h>
#include "planet.h"
#include "solver.h"
using namespace std;



int main ( int argc, char * argv[] )
{
	//NB N is now the number of initial bodies, dt is the number of integration points!
	int dim = 3, dt = atoi(argv[1]), N = atoi(argv[2]);

	double final_time =  atof(argv[3]), R0 = atof(argv[4]), M0 = 1.0;
	string particle = "particle";

	for ( int i = 0; i <= N; i++) {
		string thisparticle = particle + string(itoa(i));
	}



}