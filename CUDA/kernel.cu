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
int N = 16;

int screen_size_x = 1024;
int screen_size_y = 768;
#define BLOCK_SIZE 1024

__global__ void add_force_simple_double2(double2 *v, const double2 *r,
	const double *mass, double2 *r_r, int N, double dt)
{
	const double G = 1;

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		double _fx, _fy = 0;
		//__syncthreads();
		//#pragma unroll
		for (int j = 0; j < N; j++)
		{
			if (idx != j)
			{
				double EPS = 1E3;      // softening parameter (just to avoid infinities)
				double dx = r[j].x - r[idx].x;
				double dy = r[j].y - r[idx].y;
				double dist = sqrt(dx*dx + dy*dy);
				if (dist == 0)
				{
					dist = 0.01;
				}
				double F = (G * mass[idx] * mass[j]) / (dist*dist + EPS*EPS);

				_fx += F * dx / dist;
				_fy += F * dy / dist;
			}
			//__syncthreads();
		}
		__syncthreads();
		v[idx].x += dt * _fx / mass[idx];
		v[idx].y += dt * _fy / mass[idx];
		r_r[idx].x = v[idx].x;
		r_r[idx].y = v[idx].y;
		//r_r[idx].x += dt * v[idx].x;
		//r_r[idx].y += dt * v[idx].y;
	}
}


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

		bodies.push_back(Body(px, py, vx, vy, mass, color));
	}
	return bodies;
}

void draw_bodies(vector<Body> &bodies)
{
	for (int i = 0; i<N; i++) {
		al_draw_filled_circle((screen_size_x / 2) + (int)round(bodies[i].rx), (screen_size_y / 2) + (int)round(bodies[i].ry), bodies[i].mass / 500, bodies[i].color);
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

	// Initialise CUDA - select device
	cudaSetDevice(0);

	auto double_data_size = sizeof(double) * N;
	auto double2_data_size = sizeof(double2) * N;

	int nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
	std::cout << "size of bodies: " << N << std::endl;
	std::cout << "Blocks: " << nBlocks << ", threads per block: " << BLOCK_SIZE << std::endl;

	// Get the start time
	auto start = std::chrono::system_clock::now();
	double dt = 0.1;
	int avg_count = 10;
	for (int average_iterations = 0; average_iterations < avg_count; average_iterations++)
	{
		vector<double2> v(N);
		vector<double2> r(N);
		vector<double> mass(N);

		vector<double2> r_r(N);

		double2 *d_v;
		double2 *d_r;

		double *d_m;

		double2 *d_r_r;

		cudaMalloc((void**)&d_v, double2_data_size);
		cudaMalloc((void**)&d_r, double2_data_size);

		cudaMalloc((void**)&d_m, double_data_size);
		cudaMalloc((void**)&d_r_r, double2_data_size);

		vector<Body> bodies = startthebodies(N);

		// Get the start time
		auto current_start = std::chrono::system_clock::now();

		for (int sim_iterations = 0; sim_iterations <= 5000; sim_iterations++)
		{

			for (size_t i = 0; i < N; i++)
			{
				// Init velocity
				v[i].x = bodies[i].vx;
				v[i].y = bodies[i].vy;

				// Init positions
				r[i].x = bodies[i].rx;
				r[i].y = bodies[i].ry;

				// Init mass
				mass[i] = bodies[i].mass;
			}

			cudaMemcpyAsync(d_v, &v[0], double2_data_size, cudaMemcpyHostToDevice);
			cudaMemcpyAsync(d_r, &r[0], double2_data_size, cudaMemcpyHostToDevice);
			cudaMemcpyAsync(d_m, &mass[0], double_data_size, cudaMemcpyHostToDevice);


			add_force_simple_double2 << <nBlocks, BLOCK_SIZE >> > (d_v, d_r, d_m, d_r_r, N, dt);
			// Wait for kernel to complete
			//cudaDeviceSynchronize();
			cudaError_t error = cudaGetLastError();
			if (error != cudaSuccess)
			{
				fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(error));
				exit(-1);
			}
			// Read output buffer back to the host
			cudaMemcpyAsync(&r_r[0], d_r_r, double2_data_size, cudaMemcpyDeviceToHost);

			for (int i = 0; i < N; i++)
			{
				bodies[i].vx = r_r[i].x;
				bodies[i].vy = r_r[i].y;
				bodies[i].rx += dt * bodies[i].vx;
				bodies[i].ry += dt * bodies[i].vy;
			}

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

		std::cout << average_iterations << ", time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_total).count() << " ms" << endl;

		cudaFree(d_v);
		cudaFree(d_r);
		cudaFree(d_m);
		cudaFree(d_r_r);

	}
	// Get the end time
	auto end = std::chrono::system_clock::now();
	// Get the total time
	auto total = end - start;

	cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(total).count() / avg_count << " ms" << endl;

	al_destroy_display(display);

	return 0;
}