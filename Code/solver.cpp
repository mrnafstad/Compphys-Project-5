
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

void solver::velVerlet( int dim, int N, double final_time, bool energy)
{
	double time = 0.0;      // Sets looping variable 
	double h = final_time/(double)N;   // step length

	double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
	double acc[3];
	double acc_new[3];

	// Opening file
	FILE *fp;
	fp = fopen("VerletTest.txt", "w+");

	int counter = 0; 

	if(energy) printf("Time       Total Kinetic Energy  Total Potential Energy  Total Angular Momentum\n");


	clock_t start, finish;
	double proc_time;

	start = clock(); // starts timer


	//Starts loop
	while(time < final_time){

		fprintf(fp, "%f ", time);
		
		// Computes the gravitatonal forces between all planets in the system
		for (int j = 0; j < total_planets; j++ ) {
			planet &thisplanet = all_planets[j];

			Fx = 0; Fy = 0; Fz = 0;
			

			for (int k = 0; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx, Fy, Fz);
				}

			}

			acc[0] = Fx/thisplanet.mass; acc[1] = Fy/thisplanet.mass; acc[2] = Fz/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.position[i] += h*thisplanet.velocity[i] + 0.5*acc[i]*h*h;
			}
			
			Fx_new = 0; Fy_new = 0; Fz_new = 0;		
			
			for ( int k = 0; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx_new, Fy_new, Fz_new);
				}

			}
			acc_new[0] = Fx_new/thisplanet.mass; acc_new[1] = Fy_new/thisplanet.mass; acc_new[2] = Fz_new/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.velocity[i] += 0.5*(acc[i] + acc_new[i])*h;
			}

			fprintf(fp, "%f %f %f ", thisplanet.position[0], thisplanet.position[1], thisplanet.position[2]);
		}

		fprintf(fp, "\n");


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

	//closes file
	fclose(fp);
}


void solver::GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz){   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r));
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r));
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r));

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
