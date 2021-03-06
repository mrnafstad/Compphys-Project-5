#ifndef PLANET_H
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class planet
{
public:
	friend class solver;
	
	//Properties
	double mass;
	double position[3];
	double velocity[3];
	double potential;
	double kinetic;
	double ang_mom;
	bool bound;

	//initializers
	planet();
    planet( double M, double x, double y, double z, double vx, double vy, double vz );

	//functions
    double distance(planet otherPlanet);
    double GravitationalForce(planet otherPlanet, double Gconst);
    double Acceleration(planet otherPlanet, double Gconst);
    double KineticEnergy();
    double PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon);
};

#endif //PLANET_H