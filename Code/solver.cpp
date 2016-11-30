
#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"

solver::solver()
{
    total_planets = 0;
    radius = 100.0;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
    totalAngularMomentum = 0;

}

solver::solver( double radi )
{
    total_planets = 0;
    radius = radi;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
    totalAngularMomentum = 0;
}

void solver::add(planet newplanet)
{	
	//Adds planets to the system, with mass
    total_planets += 1;
    total_mass += newplanet.mass;
    all_planets.push_back(newplanet);
}

void solver::addM(planet newplanet)
{
	//Adds planets to the system, without mass
    total_planets +=1;
    all_planets.push_back(newplanet);
}

void solver::velVerlet( int dim, int N, double final_time, bool energy, bool stationary, bool relativity, bool MercPeri)
{
	double time = 0.0;      // Sets looping variable 
	double h = final_time/(double)N;   // step length

	planet &sun = all_planets[0];   // Fetches the sun object

	double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
	double acc[3];
	double acc_new[3];

	//Initiates helping variables for use in MercuryPerihelion()
	double rPreviousPrevious = 0;
	double rPrevious = 0;
	double previousPosition[3] = {0, 0, 0};

	// Opening files
	FILE *fp, *per;
	if(!MercPeri) fp = fopen("VerletTest.txt", "w+");
	else per = fopen("MercuryPerihelion.txt", "w+"); 


	int counter = 0; 

	if(energy) printf("Time       Total Kinetic Energy  Total Potential Energy  Total Angular Momentum\n");


	//Sets initial values for the sun if we're using center of mass as the origin
	int j, k;
	if(stationary) j = 1;
	else {
		j = 0;
		for ( int i = 1; i < total_planets; i++ ) {
			for (int k = 0; k < dim; k++ ) {
				planet &thisother = all_planets[i];
				sun.velocity[k] -= thisother.velocity[k]*thisother.mass/sun.mass;
			}
		}
	}

	clock_t start, finish;
	double proc_time;

	start = clock(); // starts timer


	//Starts loop
	while(time < final_time){

		//Writes to file if MercPeri = false
		if(!MercPeri) fprintf(fp, "%f ", time);

		if(stationary) j = 1;
		else j = 0;
		
		// Computes the gravitatonal forces between all planets in the system
		for ( j; j < total_planets; j++ ) {
			planet &thisplanet = all_planets[j];

			// Computers the perihelion of the orbit of Mercury
			if(MercPeri) MercuryPerihelion(thisplanet, sun, rPreviousPrevious, rPrevious, previousPosition, time, per);

			Fx = 0; Fy = 0; Fz = 0;
			

			if (stationary) {
				GravitationalForce(thisplanet, sun, Fx, Fy, Fz, relativity);
				k = 1;
			}
			else k = 0;

			for ( k; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx, Fy, Fz, relativity );
				}

			}

			acc[0] = Fx/thisplanet.mass; acc[1] = Fy/thisplanet.mass; acc[2] = Fz/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.position[i] += h*thisplanet.velocity[i] + 0.5*acc[i]*h*h;
			}
			
			Fx_new = 0; Fy_new = 0; Fz_new = 0;
			if (stationary) {
				GravitationalForce(thisplanet, sun, Fx_new, Fy_new, Fz_new, relativity);
				k = 1;
			}
			else k = 0;			
			
			for ( k; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx_new, Fy_new, Fz_new, relativity );
				}

			}
			acc_new[0] = Fx_new/thisplanet.mass; acc_new[1] = Fy_new/thisplanet.mass; acc_new[2] = Fz_new/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.velocity[i] += 0.5*(acc[i] + acc_new[i])*h;
			}

			// Prints position to file if not computing Mercury Perihelion
			if(!MercPeri) fprintf(fp, "%f %f %f ", thisplanet.position[0], thisplanet.position[1], thisplanet.position[2]);
		}

		if(!MercPeri) fprintf(fp, "\n");


		// Prints energies to screen if required
		if(energy){
			if(counter == 0){
				KineticEnergySystem();
				PotentialEnergySystem(0.0);
				AngularMomentumSystem();
				printf("%f      %e        %e            %e\n", time, totalKinetic, totalPotential, totalAngularMomentum);
			}
		}

		//Makes sure that the energy above only prints every 1/10th iteration 
		counter += 1;
		if(N/counter == 10) counter = 0;
		time += h;

	}
	finish = clock();  // stopping timer
	proc_time = ( (double) (finish - start)/CLOCKS_PER_SEC);
	printf("Time spent on algorithm: %f seconds\n", proc_time);

	//closes files
	if(!MercPeri) fclose(fp);
	else fclose(per);
}


void solver::ForwardEuler( int dim, int N, double final_time, bool relativity)
{
	double time = 0.0;     // Sets the looping variable
	double h = final_time/(double)N;  // step length

	//Fetches the sun and earth objects
	planet &earth = all_planets[1];   
	planet &sun = all_planets[0];

	double Fx, Fy, Fz;
	double acc[3];

	// Opens file to write positions to
	FILE *fp;
	fp = fopen("EulerTest.txt", "w+");
	fprintf(fp, "%f %f %f %f\n", time, earth.position[0], earth.position[1], earth.position[2]);

	clock_t start, finish;
	double proc_time;

	start = clock(); // Starting timer
	while(time < final_time){

		//Computing forces and acceleration
		Fx = 0, Fy = 0, Fz = 0;
		GravitationalForce(earth, sun, Fx, Fy, Fz, relativity);

		acc[0] = Fx/earth.mass; acc[1] = Fy/earth.mass; acc[2] = Fz/earth.mass;

		//Forward Euler
		for(int i = 0; i < dim; i++){

			earth.position[i] += h*earth.velocity[i];
			earth.velocity[i] += h*acc[i];
		}
		fprintf(fp, "%f %f %f %f\n", time, earth.position[0], earth.position[1], earth.position[2]);
			
		time += h; // updates loop variable
	}

	finish = clock(); // stops timer
	proc_time = ((double) (finish - start)/CLOCKS_PER_SEC);

	printf("Time spent on algorithm: %f seconds.\n", proc_time);

	fclose(fp);
		
}

void solver::GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, bool relativity){   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r));
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r));
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r));

    if(relativity){
    	//Adds the relativistic correction if required
    	double angmom2 = current.AngularMomentum()*current.AngularMomentum()/(current.mass*current.mass);
    	double c = 63197.8;  //Speed of light AU/yr
    	Fx *= (1 + 3*angmom2/(r*r*c*c));
    	Fy *= (1 + 3*angmom2/(r*r*c*c));
    	Fz *= (1 + 3*angmom2/(r*r*c*c));
    }
}

void solver::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.kinetic = Current.KineticEnergy();
        totalKinetic += Current.kinetic;
    }
}

void solver::PotentialEnergySystem(double epsilon)
{
    totalPotential = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.potential = 0;
    }
    for(int nr1=0;nr1<total_planets;nr1++){
        planet &Current = all_planets[nr1];
        for(int nr2=nr1+1;nr2<total_planets;nr2++){
            planet &Other = all_planets[nr2];
            Current.potential += Current.PotentialEnergy(Other,G,epsilon);
            Other.potential += Other.PotentialEnergy(Current,G,epsilon);
        	totalPotential += Current.potential;
        }
    }
}

void solver::AngularMomentumSystem(){
	totalAngularMomentum = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.ang_mom = Current.AngularMomentum();
        totalAngularMomentum += Current.ang_mom;
    }
}

void solver::MercuryPerihelion(planet &thisplanet, planet &sun, double &rPreviousPrevious, double &rPrevious, double previousPosition[], double time, FILE *per){
	double rCurrent = thisplanet.distance(sun);

	//Checks if the planet is at perihelion by checking if the previous distance to the sun
	// is shorter than the ones in the previousprevious and the current iterations
	if( rCurrent > rPrevious && rPrevious < rPreviousPrevious){
		double x = previousPosition[0];
		double y = previousPosition[1];
		printf("Time: %f, Perihelion angle: %f rad = %f ''\n", time, atan2(y,x), atan2(y,x)*648000/M_PI);
		fprintf(per, "Time: %f, Perihelion angle: %f rad = %f ''\n", time, atan2(y,x), atan2(y,x)*648000/M_PI);
	}

	//Updates variables 
	rPreviousPrevious = rPrevious;
	rPrevious = rCurrent;
	for(int i = 0; i < 3; i++){
		previousPosition[i] = thisplanet.position[i];
	}
}