#pragma once
using namespace std;
#include <math.h>
#include <allegro5/allegro_color.h>
#include <iostream>
class Body 
{
private:
	 const double G = 1;   // gravitational constant
	 const double solarmass = 1.98892e30;
public:
	double rx, ry;       // holds the cartesian positions
	double vx, vy;       // velocity components 
	double fx, fy;       // force components
	double mass;         // mass
	ALLEGRO_COLOR color; // color of body

	Body();

	// create and initialize a new Body
	Body(double rx, double ry, double vx, double vy, double mass, ALLEGRO_COLOR color);

	// update the velocity and position using a timestep dt
	void update(double dt);

	// returns the distance between two bodies
	double distanceTo(Body &b);

	// set the force to 0 for the next iteration
	void resetForce();

	// compute the net force acting between the body a and b, and
	// add to the net force acting on a
	void addForce(Body &b);
};