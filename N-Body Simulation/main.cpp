#include "body.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <allegro5/allegro_color.h>
#include <iostream>
using namespace std;

int N = 250;
vector<Body*> bodies; //Body bodies[1000];

double random()
{
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

//double exp(double lambda) 
//{
//	return -log(1 - random()) / lambda;
//}

int sgn(double d) {
	return d<-DBL_EPSILON ? -1 : d>DBL_EPSILON;
}

//the bodies are initialized in circular orbits around the central mass.
//This is just some physics to do that
double circlev(double rx, double ry)
{
	double solarmass = 1.98892e30;
	double r2 = sqrt(rx*rx + ry*ry);
	double numerator = (6.67e-11)*1e6*solarmass;
	return sqrt(numerator / r2);
}

//Initialize N bodies with random positions and circular velocities
void startthebodies(int N)
{
	double radius = 1e18;        // radius of universe
	double solarmass = 1.98892e30;
	for (int i = 0; i < N; i++) {
		double px = 1e18*exp(-1.8)*(.5 - random());
		double py = 1e18*exp(-1.8)*(.5 - random());
		double magv = circlev(px, py);

		double absangle = atan(abs(py / px));
		double thetav = M_PI / 2 - absangle;
		double phiv = random() * M_PI;
		double vx = -1 * sgn(py)*cos(thetav)*magv;
		double vy = sgn(px)*sin(thetav)*magv;
		// Orient a random 2D circular orbit
		if (random() <= .5) {
			vx = -vx;
			vy = -vy;
		}

		double mass = random() * solarmass * 10 + 1e20;
		// Color the masses in green gradients by mass
		int red = (int)floor(mass * 254 / (solarmass * 10 + 1e20));
		int blue = (int)floor(mass * 254 / (solarmass * 10 + 1e20));
		int green = 255;
		ALLEGRO_COLOR color = al_map_rgb(red, green, blue);
		// put a heavy body in the center
		if (i == 0)
			bodies.push_back(new Body(0, 0, 0, 0, 1e6*solarmass, color));

		bodies.push_back(new Body(px, py, vx, vy, mass, color));
	}
}

//Use the method in Body to reset the forces, then add all the new forces
void addforces(int N)
{
	for (int i = 0; i < N; i++) {
		bodies[i]->resetForce();
		//Notice-2 loops-->N^2 complexity
		for (int j = 0; j < N; j++) {
			if (i != j) bodies[i]->addForce(*bodies[j]);
		}
	}
	//Then, loop again and update the bodies using timestep dt
	for (int i = 0; i < N; i++) {
		bodies[i]->update(1e11);
	}
}

//void init_window()
//{
//	al_init();
//	al_create_display(1024, 768);
//}

void draw_bodies()
{
	//std::cout << random() << std::endl;
	//std::cout << bodies[0]->rx << " " << bodies[0]->ry << std::endl;
	for (int i = 0; i<N; i++) { 
		al_draw_circle((1024/2) + (int)round(bodies[i]->rx / 1e18), (768 / 2) + (int)round(bodies[i]->ry / 1e18),1.0f, bodies[i]->color, 0.75f);
		//g.fillOval((int)round(bodies[i].rx * 250 / 1e18), (int)round(bodies[i].ry * 250 / 1e18), 8, 8);
	}
}

int main(int argc, char **argv)
{

	ALLEGRO_DISPLAY *display = NULL;

	if (!al_init()) {
		fprintf(stderr, "failed to initialize allegro!\n");
		return -1;
	}

	display = al_create_display(1024, 768);
	if (!display) {
		fprintf(stderr, "failed to create display!\n");
		return -1;
	}

	al_clear_to_color(al_map_rgb(0, 0, 0));

	al_flip_display();
	al_init_primitives_addon();

	//init_window();
	startthebodies(N);

	while (true)
	{
		addforces(N);
		draw_bodies();
		al_flip_display();
		al_clear_to_color(al_map_rgb(0, 0, 0));
	}

	al_destroy_display(display);

	return 0;
}