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
    void ForwardEuler(int dim, int N, double final_time, bool relativity);
    void velVerlet(int dim, int N, double final_time, bool energy, bool stationary, bool relativity, bool MercPeri);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, bool relativity);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    void AngularMomentumSystem();
    void MercuryPerihelion(planet &thisplanet, planet &sun, double &rPreviousPrevious, double &rPrevious, double previousPosition[], double time, FILE *per);
};

#endif //SOLVER_H