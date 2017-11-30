// Translated into C++ from Java, based on the code available at: http://physics.princeton.edu/~fpretori/Nbody/

#include "body.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <allegro5/allegro_color.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_ttf.h>
#include <iostream>
#include <chrono>
using namespace std;

int N = 1024;
vector<Body> bodies; //Body bodies[1000];

int screen_size_x = 1024;
int screen_size_y = 768;

double random()
{
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
int sgn(double d) {
	return d<-DBL_EPSILON ? -1 : d>DBL_EPSILON;
}

// The bodies are initialized in circular orbits around the central mass.
// This is just some physics to do that
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
	//double solarmass = 1.98892e30;
	for (int i = 0; i < N; i++) {
		double px = (rand() % screen_size_x) - screen_size_x / 2; //exp(-1.8)*(.5 - random());
		double py = (rand() % screen_size_y) - screen_size_y / 2; // exp(-1.8)*(.5 - random());
		double magv = circlev(px, py);

		double absangle = atan(abs(py / px));
		double thetav = M_PI / 2 - absangle;
		double phiv = random() * M_PI;
		double vx = -1 * sgn(py)*cos(thetav)*magv;
		double vy = sgn(px)*sin(thetav)*magv;

		vx = 0;
		vy = 0;

		// Orient a random 2D circular orbit
		//if (random() <= .5) {
		//	vx = -vx;
		//	vy = -vy;
		//}

		double mass = rand() % 20 + 1;  //random() * solarmass * 10 + 1e20;
		// Color the masses in green gradients by mass
		int red = (int)floor(255);
		int blue = (int)floor(255);
		int green = 255;
		ALLEGRO_COLOR color = al_map_rgb(red, green, blue);
		// put a heavy body in the center
		if (i == 0)
			bodies.push_back(Body(0, 0, 0, 0, 100, color));

		bodies.push_back(Body(px, py, vx, vy, mass, color));
	}
}

//Use the method in Body to reset the forces, then add all the new forces
void addforces(int N)
{
	for (int i = 0; i < N; i++) {
		bodies[i].resetForce();
		//Notice-2 loops-->N^2 complexity
		for (int j = 0; j < N; j++) {
			if (i != j) bodies[i].addForce(bodies[j]);
		}
		bodies[i].update(1e11);
	}
	//Then, loop again and update the bodies using timestep dt
	//for (int i = 0; i < N; i++) {
	//	bodies[i]->update(1e11);
	//}
}

void draw_bodies()
{
	for (int i = 0; i<N; i++) { 
		al_draw_circle((screen_size_x /2) + (int)round(bodies[i].rx), (screen_size_y / 2) + (int)round(bodies[i].ry),1.0f, bodies[i].color, 0.75f);
	}
}

int main(int argc, char **argv)
{

	ALLEGRO_DISPLAY *display = NULL;

	if (!al_init()) {
		fprintf(stderr, "failed to initialize allegro!\n");
		return -1;
	}

	display = al_create_display(screen_size_x, screen_size_y);
	if (!display) {
		fprintf(stderr, "failed to create display!\n");
		return -1;
	}

	al_clear_to_color(al_map_rgb(0, 0, 0));

	al_flip_display();
	al_init_primitives_addon();
	al_init_font_addon(); // initialize the font addon
	al_init_ttf_addon();// initialize the ttf (True Type Font) addon

	ALLEGRO_FONT *font = al_load_ttf_font("../Consolas.ttf", 24, 0);

	if (!font) {
		fprintf(stderr, "Could not load 'Consolas.ttf'.\n");
		return -1;
	}



	// Get the start time
	auto start = std::chrono::system_clock::now();
	int avg_count = 100;

	for (int average_iterations = 0; average_iterations <= avg_count; average_iterations++)
	{
		bodies.clear();
		startthebodies(N);
		// Get the start time
		auto current_start = std::chrono::system_clock::now();

		for (int sim_iterations = 0; sim_iterations <= 100000; sim_iterations++)
		{
			addforces(N);
			draw_bodies();
			al_flip_display();
			al_clear_to_color(al_map_rgb(0, 0, 0));

			char buffer[30];
			_itoa_s(sim_iterations, buffer, 10);
			string s = "Simulation iterations: ";
			s += buffer;
			al_draw_text(font, al_map_rgb(255, 255, 255), 0, 0, ALLEGRO_ALIGN_LEFT, s.c_str());
		}
		// Get the end time
		auto current_end = std::chrono::system_clock::now();
		// Get the total time
		auto current_total = current_end - current_start;

		cout << average_iterations << ", time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_total).count() << " ms" << endl;
	}
	// Get the end time
	auto end = std::chrono::system_clock::now();
	// Get the total time
	auto total = end - start;

	cout << "Total time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(total).count() / avg_count << " ms" << endl;

	al_destroy_display(display);

	return 0;
}