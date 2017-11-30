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
#include <thread>
using namespace std;

int N = 2048;
//vector<Body> bodies; //Body bodies[1000];
int num_threads;
int screen_size_x = 1024;
int screen_size_y = 768;

//Initialize N bodies with random positions and circular velocities
vector<Body> startthebodies(int N)
{
	vector<Body> bodies;
	for (int i = 0; i < N; i++) {
		// Initialise position
		double px = (rand() % screen_size_x) - screen_size_x / 2;
		double py = (rand() % screen_size_y) - screen_size_y / 2;

		// Initialise velocity
		double vx = 0;
		double vy = 0;

		double mass = rand() % 1000 + 1;
		// Color the masses in blue gradients by mass
		int red = (int)floor(mass * 254);
		int blue = 255;
		int green = (int)floor(mass * 254);
		ALLEGRO_COLOR color = al_map_rgb(red, green, blue);

		if (i == 0)
			std::cout << "x: " << px << ", y:" << py << std::endl;

		bodies.push_back(Body(px, py, vx, vy, mass, color));
	}
	return bodies;
}

//Use the method in Body to reset the forces, then add all the new forces
void addforces(vector<Body> &bodies, int N)
{
	#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for (int i = 0; i < N; i++) {
		bodies[i].resetForce();
		//Notice-2 loops-->N^2 complexity
		for (int j = 0; j < N; j++) {
			if (i != j) bodies[i].addForce(bodies[j]);
		}
	}
	//#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	//Then, loop again and update the bodies using timestep dt
	for (int i = 0; i < N; i++) {
		bodies[i].update(1);
	}
}

void draw_bodies(vector<Body> &bodies)
{
	for (int i = 0; i<N; i++) { 
		al_draw_filled_circle((screen_size_x / 2) + (int)round(bodies[i].rx), (screen_size_y / 2) + (int)round(bodies[i].ry), bodies[i].mass / 1000, bodies[i].color);
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
	num_threads = thread::hardware_concurrency();
	// Get the start time
	auto start = std::chrono::system_clock::now();
	int avg_count = 100;

	for (int average_iterations = 0; average_iterations <= avg_count; average_iterations++)
	{
		//bodies.clear();
		vector<Body> bodies = startthebodies(N);
		// Get the start time
		auto current_start = std::chrono::system_clock::now();
		std::cout << "size of bodies: " << bodies.size() << std::endl;
		for (int sim_iterations = 0; sim_iterations < 1000; sim_iterations++)
		{
			addforces(bodies, N);
			draw_bodies(bodies);
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