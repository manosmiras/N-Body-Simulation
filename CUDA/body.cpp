#include "body.h"

Body::Body()
{
}

Body::Body(double rx, double ry, double vx, double vy, double mass, ALLEGRO_COLOR color)
{
	this->rx = rx;
	this->ry = ry;
	this->vx = vx;
	this->vy = vy;
	this->mass = mass;
	this->color = color;
	this->fx = 0;
	this->fy = 0;
}

// update the velocity and position using a timestep dt
void Body::update(double dt) 
{
	vx += dt * fx / mass;
	vy += dt * fy / mass;
	rx += dt * vx;
	ry += dt * vy;
}

// returns the distance between two bodies
double Body::distanceTo(Body &b) 
{
	double dx = rx - b.rx;
	double dy = ry - b.ry;
	return sqrt(dx*dx + dy*dy);
}

// set the force to 0 for the next iteration
void Body::resetForce() 
{
	fx = 0.0;
	fy = 0.0;
}

// compute the net force acting between the body a and b, and
// add to the net force acting on a
void Body::addForce(Body &b) 
{
	//Body a = *this;
	double EPS = 1E3;      // softening parameter (just to avoid infinities)
	double dx = b.rx - this->rx;
	double dy = b.ry - this->ry;
	double dist = sqrt(dx*dx + dy*dy);
	if (dist == 0)
	{
		dist = 0.01;
	}
	//	std::cout << "dist is 0" << std::endl;
	double F = (G * this->mass * b.mass) / (dist*dist + EPS * EPS);
	//std::cout << F << std::endl;
	this->fx += F * dx / dist;
	this->fy += F * dy / dist;
}
