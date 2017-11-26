// Translated into C++ from Java, based on the code available at: http://physics.princeton.edu/~fpretori/Nbody/
#define _USE_MATH_DEFINES
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "body.h"
#include <vector>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <allegro5/allegro_color.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_ttf.h>
#include <iostream>
#include <chrono>
#include <string>
using namespace std;
const double G = 6.673e-11;   // gravitational constant
int N = 1000;
vector<Body*> bodies; //Body bodies[1000];

int screen_size_x = 1024;
int screen_size_y = 768;
#define BLOCK_SIZE 1024

struct force
{
	double fx;
	double fy;
};

__global__ void add_force(Body *bodies, force *f, const double G, int N, double dt)
{
	// Get block index
	unsigned int block_idx = blockIdx.x;
	// Get thread index
	unsigned int thread_idx = threadIdx.x;
	// Get the number of threads per block
	unsigned int block_dim = blockDim.x;
	// Get the thread's unique ID = (block_idx * block_dim) + thread_idx;
	unsigned int idx = (block_idx * block_dim) + thread_idx;

	// Reset forces
	bodies[idx].fx = 0.0;
	bodies[idx].fy = 0.0;

	for (int j = 0; j < N; j++) 
	{
		if (idx != j)
		{
			double EPS = 3E4;      // softening parameter (just to avoid infinities)
			double dx = bodies[j].rx - bodies[idx].rx;
			double dy = bodies[j].ry - bodies[idx].ry;
			double dist = sqrt(dx*dx + dy*dy);
			double F = (G * bodies[idx].mass * bodies[j].mass) / (dist*dist + EPS*EPS);

			//fx[idx] = bodies[idx].fx;
			//fy[idx] = bodies[j].fx;
			f[idx].fx += F * dx / dist;
			f[idx].fy += F * dy / dist;
			//bodies_return[idx].fx += F * dx / dist;
			//bodies_return[idx].fy += F * dy / dist;
		}
	}

	// Update velocity
	//bodies_return[idx].vx += dt * bodies_return[idx].fx / bodies_return[idx].mass;
	//bodies_return[idx].vy += dt * bodies_return[idx].fy / bodies_return[idx].mass;
	//bodies_return[idx].rx += dt * bodies_return[idx].vx;
	//bodies_return[idx].ry += dt * bodies_return[idx].vy;
	//c->fx += F * dx / dist;
	//c->fy += F * dy / dist;
}


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

// Initialize N bodies with random positions and circular velocities
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

// Use the method in Body to reset the forces, then add all the new forces
void addforces(int N)
{
	for (int i = 0; i < N; i++) {
		bodies[i]->resetForce();
		// Notice-2 loops-->N^2 complexity
		for (int j = 0; j < N; j++) {
			if (i != j) bodies[i]->addForce(*bodies[j]);
		}
	}
	// Then, loop again and update the bodies using timestep dt
	for (int i = 0; i < N; i++) {
		bodies[i]->update(1e11);
	}
}

void draw_bodies()
{
	for (int i = 0; i<N; i++) {
		al_draw_circle((screen_size_x / 2) + (int)round(bodies[i]->rx / 1e18), (screen_size_y / 2) + (int)round(bodies[i]->ry / 1e18), 1.0f, bodies[i]->color, 0.75f);
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
	startthebodies(N);

	// Initialise CUDA - select device
	cudaSetDevice(0);

	// Create host memory
	auto body_data_size = sizeof(Body) * N;
	auto double_data_size = sizeof(force) * N;
	vector<force> f(N); // Output array
	//vector<double> fy(N); // Output array

	// Declare buffers
	Body *buffer_A;
	force *buffer_B;//, *buffer_C;

	// Initialise buffers
	cudaMalloc((void**)&buffer_A, body_data_size);
	cudaMalloc((void**)&buffer_B, double_data_size);
	//cudaMalloc((void**)&buffer_C, double_data_size);

	int nBlocks = N / BLOCK_SIZE;

	// Get the start time
	auto start = std::chrono::system_clock::now();

	for (int sim_iterations = 0; sim_iterations <= 5000; sim_iterations++)
	{
		//addforces(N);
		cudaMemcpy(buffer_A, &bodies[0], body_data_size, cudaMemcpyHostToDevice);

		// Run kernel with one thread for each element
		add_force <<<nBlocks, BLOCK_SIZE>>>(buffer_A, buffer_B, G, N, 1e11);
		//kernel << <1, 1 >> > ();
		// Wait for kernel to complete
		cudaDeviceSynchronize();

		// Read output buffers back to the host
		cudaMemcpy(&f[0], buffer_B, double_data_size, cudaMemcpyDeviceToHost);
		//cudaMemcpy(&fy[0], buffer_C, double_data_size, cudaMemcpyDeviceToHost);

		for (int i = 0; i < N; i++) 
		{
			bodies[i]->fx = f[i].fx;
			bodies[i]->fy = f[i].fy;

			bodies[i]->update(1e11);
		}

		draw_bodies();
		al_flip_display();
		al_clear_to_color(al_map_rgb(0, 0, 0));
		char buffer[30];
		itoa(sim_iterations, buffer, 10);
		string s = "Simulation iterations: ";
		s += buffer;
		al_draw_text(font, al_map_rgb(255, 255, 255), 0, 0, ALLEGRO_ALIGN_LEFT, s.c_str());
	}

	// Get the end time
	auto end = std::chrono::system_clock::now();
	// Get the total time
	auto total = end - start;

	cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(total).count() << " ms" << endl;

	al_destroy_display(display);

	// Clean up resources
	cudaFree(buffer_A);
	cudaFree(buffer_B);
	//cudaFree(buffer_C);

	return 0;
}