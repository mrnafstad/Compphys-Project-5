#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>

using std::vector;

class solver
{
public:
	friend class planet;

	//properties
	double radius, G;
	vector<planet> all_planets;
	int total_planets;
	double total_mass;
	double totalKinetic;
	double totalPotential;
	double totalAngularMomentum;


	//initializer
	solver();
	solver( double radii );

	//functions;
	void Gravitiationalconstant();
	void add(planet newPlanet);
	void addM(planet newPlanet);
    void velVerlet(int dim, double dt, double final_time, int N, bool energy);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    void AngularMomentumSystem();
};

#endif //SOLVER_H