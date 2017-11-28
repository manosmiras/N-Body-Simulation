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
int N = 1024;
vector<Body*> bodies; //Body bodies[1000];

int screen_size_x = 1024;
int screen_size_y = 768;
#define BLOCK_SIZE 512

//struct force
//{
//	double fx;
//	double fy;
//};

struct Body_struct
{
	double vx;
	double vy;
	double rx;
	double ry;
	double mass;
};

__global__ void add_force_block_doubles(double *vx, double *vy, const double *rx, const double *ry,
	const double *mass, double *r_rx, double *r_ry, const double G, int N, double dt)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		// Reset forces
		double fx, fy = 0.0;

		for (int tile = 0; tile < gridDim.x; tile++)
		{
			__shared__ double s_rx[BLOCK_SIZE];
			__shared__ double s_ry[BLOCK_SIZE];
			double t_rx = rx[tile * blockDim.x + threadIdx.x];
			double t_ry = ry[tile * blockDim.x + threadIdx.x];
			s_rx[threadIdx.x] = t_rx;
			s_ry[threadIdx.x] = t_ry;
			__syncthreads();
			for (int j = 0; j < BLOCK_SIZE; j++)
			{
				if (idx != j)
				{
					double EPS = 3E4;      // softening parameter (just to avoid infinities)
					double dx = s_rx[j] - rx[idx];
					double dy = s_ry[j] - ry[idx];
					double dist = sqrt(dx*dx + dy*dy);
					double F = (G * mass[idx] * mass[j]) / (dist*dist + EPS*EPS);

					fx += F * dx / dist;
					fy += F * dy / dist;
				}
				__syncthreads();
			}
			//__syncthreads();
		}

		// Calculate velocity and integrate
		vx[idx] += dt * fx / mass[idx];
		vy[idx] += dt * fy / mass[idx];
		r_rx[idx] += dt * vx[idx];
		r_ry[idx] += dt * vy[idx];
	}
}

__global__ void add_force(const Body *bodies, double *vx, double *vy, const double G, int N, double dt)
{
	// Get block index
	unsigned int block_idx = blockIdx.x;
	// Get thread index
	unsigned int thread_idx = threadIdx.x;
	// Get the number of threads per block
	unsigned int block_dim = blockDim.x;
	// Get the thread's unique ID = (block_idx * block_dim) + thread_idx;
	unsigned int idx = (block_idx * block_dim) + thread_idx;
	if (idx < N)
	{
		// Reset forces
		double fx, fy = 0.0;
		//bodies[idx].fx = 0.0;
		//bodies[idx].fy = 0.0;

		for (int tile = 0; tile < gridDim.x; tile++)
		{
			__shared__ double s_rx[BLOCK_SIZE];
			__shared__ double s_ry[BLOCK_SIZE];
			double t_rx = bodies[tile * blockDim.x + threadIdx.x].rx;
			double t_ry = bodies[tile * blockDim.x + threadIdx.x].ry;
			s_rx[threadIdx.x] = t_rx;
			s_ry[threadIdx.x] = t_ry;
			__syncthreads();
			for (int j = 0; j < BLOCK_SIZE; j++)
			{
				if (idx != j)
				{
					double EPS = 3E4;      // softening parameter (just to avoid infinities)
					double dx = s_rx[j] - bodies[idx].rx;
					double dy = s_ry[j] - bodies[idx].ry;
					double dist = sqrt(dx*dx + dy*dy);
					double F = (G * bodies[idx].mass * bodies[j].mass) / (dist*dist + EPS*EPS);

					// is this right? probably
					//f[idx].fx += F * dx / dist;
					//f[idx].fy += F * dy / dist;

					fx += F * dx / dist;
					fy += F * dy / dist;
				}
				
			}
			__syncthreads();
		}
		
		// Calculate velocity
		vx[idx] += dt * fx / bodies[idx].mass;
		vy[idx] += dt * fy / bodies[idx].mass;
	}
	//for (int j = 0; j < N; j++) 
	//{
	//	if (idx != j)
	//	{
	//		double EPS = 3E4;      // softening parameter (just to avoid infinities)
	//		double dx = bodies[j].rx - bodies[idx].rx;
	//		double dy = bodies[j].ry - bodies[idx].ry;
	//		double dist = sqrt(dx*dx + dy*dy);
	//		double F = (G * bodies[idx].mass * bodies[j].mass) / (dist*dist + EPS*EPS);

	//		//fx[idx] = bodies[idx].fx;
	//		//fy[idx] = bodies[j].fx;
	//		f[idx].fx += F * dx / dist;
	//		f[idx].fy += F * dy / dist;
	//		//bodies_return[idx].fx += F * dx / dist;
	//		//bodies_return[idx].fy += F * dy / dist;
	//	}
	//}

	// Update velocity
	//bodies_return[idx].vx += dt * bodies_return[idx].fx / bodies_return[idx].mass;
	//bodies_return[idx].vy += dt * bodies_return[idx].fy / bodies_return[idx].mass;
	//bodies_return[idx].rx += dt * bodies_return[idx].vx;
	//bodies_return[idx].ry += dt * bodies_return[idx].vy;
	//c->fx += F * dx / dist;
	//c->fy += F * dy / dist;
}

__global__ void add_force_simple(const Body *bodies, Body *returned_bodies, const double G, int N, double dt)//(Body *bodies, double *vx, double *vy, const double G, int N, double dt)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		double fx, fy = 0;
		returned_bodies[idx].fx = 0;
		returned_bodies[idx].fy = 0;
		for (int j = 0; j < N; j++)
		{
			if (idx != j)
			{
				double EPS = 3E4;      // softening parameter (just to avoid infinities)
				double dx = bodies[j].rx - bodies[idx].rx;
				double dy = bodies[j].ry - bodies[idx].ry;
				double dist = sqrt(dx*dx + dy*dy);
				double F = (G * bodies[idx].mass * bodies[j].mass) / (dist*dist + EPS*EPS);

				fx += F * dx / dist;
				fy += F * dy / dist;
				//returned_bodies[idx].fx += fx;
				//returned_bodies[idx].fy += fy;
				//fx[idx] = bodies[idx].fx;
				//fy[idx] = bodies[j].fx;
				//f[idx].fx += F * dx / dist;
				//f[idx].fy += F * dy / dist;
				//bodies_return[idx].fx += F * dx / dist;
				//bodies_return[idx].fy += F * dy / dist;
			}
		}
		//bodies[idx].vx += dt * fx / bodies[idx].mass;
		//bodies[idx].vy += dt * fy / bodies[idx].mass;
		//bodies[idx].rx += dt * bodies[idx].vx;
		//bodies[idx].ry += dt * bodies[idx].vy;

		// Update velocity and integrate
		//returned_bodies[idx].vx += dt * returned_bodies[idx].fx / returned_bodies[idx].mass;
		//returned_bodies[idx].vy += dt * returned_bodies[idx].fy / returned_bodies[idx].mass;
		//returned_bodies[idx].rx += dt * returned_bodies[idx].vx;
		//returned_bodies[idx].ry += dt * returned_bodies[idx].vy;
		//vx[idx] += dt * fx / bodies[idx].mass;
		//vy[idx] += dt * fy / bodies[idx].mass;
	}

	//bodies_return[idx].vx += dt * bodies_return[idx].fx / bodies_return[idx].mass;
	//bodies_return[idx].vy += dt * bodies_return[idx].fy / bodies_return[idx].mass;
	//bodies_return[idx].rx += dt * bodies_return[idx].vx;
	//bodies_return[idx].ry += dt * bodies_return[idx].vy;
	//c->fx += F * dx / dist;
	//c->fy += F * dy / dist;
}

__global__ void add_force_simple_struct(Body_struct *bodies, double *r_rx, double *r_ry, const double G, int N, double dt)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		double _fx, _fy = 0;
		for (int j = 0; j < N; j++)
		{
			if (idx != j)
			{
				double EPS = 3E4;      // softening parameter (just to avoid infinities)
				double dx = bodies[j].rx - bodies[idx].rx;
				double dy = bodies[j].ry - bodies[idx].ry;
				double dist = sqrt(dx*dx + dy*dy);
				double F = (G * bodies[idx].mass * bodies[j].mass) / (dist*dist + EPS*EPS);

				_fx += F * dx / dist;
				_fy += F * dy / dist;
			}
			__syncthreads();
		}
		bodies[idx].vx += dt * _fx / bodies[idx].mass;
		bodies[idx].vy += dt * _fy / bodies[idx].mass;
		r_rx[idx] += dt * bodies[idx].vx;
		r_ry[idx] += dt * bodies[idx].vy;
	}
}



__global__ void add_force_simple_doubles(double *vx, double *vy, const double *rx, const double *ry,
	const double *mass, double *r_rx, double *r_ry, const double G, int N, double dt)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < N)
	{
		double _fx, _fy = 0;
		//__syncthreads();
		#pragma unroll
		for (int j = 0; j < N; j++)
		{
			if (idx != j)
			{
				double EPS = 3E4;      // softening parameter (just to avoid infinities)
				double dx = rx[j] - rx[idx];
				double dy = ry[j] - ry[idx];
				double dist = sqrt(dx*dx + dy*dy);
				double F = (G * mass[idx] * mass[j]) / (dist*dist + EPS*EPS);

				_fx += F * dx / dist;
				_fy += F * dy / dist;
				//
			}
			//__syncthreads();
		}
		__syncthreads();
		vx[idx] += dt * _fx / mass[idx];
		vy[idx] += dt * _fy / mass[idx];
		r_rx[idx] += dt * vx[idx];
		r_ry[idx] += dt * vy[idx];
	}
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
		//al_draw_circle((screen_size_x / 2) + (int)round(rx[i] / 1e18), (screen_size_y / 2) + (int)round(ry[i] / 1e18), 1.0f, bodies[i]->color, 0.75f);
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
	//auto body_data_size = sizeof(Body*) * N;
	//auto returned_body_data_Size = sizeof(Body) * N;
	//auto force_data_size = sizeof(force*) * N;

	auto double_data_size = sizeof(double) * N;
	//auto body_struct_data_size = sizeof(Body_struct) * N;

	//vector<Body> returned_bodies(N);
	
	//for (size_t i = 0; i < N; i++)
	//{
	//	returned_bodies.push_back(*bodies[i]);
	//}

	vector<double> vx(N);
	vector<double> vy(N);
	vector<double> rx(N);
	vector<double> ry(N);
	vector<double> fx(N);
	vector<double> fy(N);
	vector<double> mass(N);
	//
	for (size_t i = 0; i < N; i++)
	{
		// Init velocity
		vx[i] = bodies[i]->vx;
		vy[i] = bodies[i]->vy;
		// Init positions
		rx[i] = bodies[i]->rx;
		ry[i] = bodies[i]->ry;
		// Init forces
		//fx[i] = bodies[i]->fx;
		//fy[i] = bodies[i]->fy;
		// Init mass
		mass[i] = bodies[i]->mass;
	}
	//vector<Body_struct> bs(N);
	//for (size_t i = 0; i < N; i++)
	//{
	//	// Init velocity
	//	bs[i].vx = bodies[i]->vx;
	//	bs[i].vy = bodies[i]->vy;
	//	// Init positions
	//	bs[i].rx = bodies[i]->rx;
	//	bs[i].ry = bodies[i]->ry;
	//	// Init mass
	//	bs[i].mass = bodies[i]->mass;
	//}

	//vector<force*> f(N); // Output array
	vector<double> r_rx(N);
	vector<double> r_ry(N);
	//
	//for (size_t i = 0; i < N; i++)
	//{
	//	f[i] = (force*)malloc(force_data_size);
	//	f[i]->fx = 0.0;
	//	f[i]->fy = 0.0;
	//}
	
	// Declare buffers
	//Body *d_a;// = (Body*)malloc(body_data_size);
	////force *d_b;
	////double *d_b;
	//Body *d_b;
	//double *d_c;

	double *d_vx;
	double *d_vy;
	double *d_rx;
	double *d_ry;
	double *d_m;
	double *d_r_rx;
	double *d_r_ry;
	//Body_struct *d_bs;

	// Initialise buffers
	//cudaMalloc((void**)&d_a, body_data_size);
	//cudaMalloc((void**)&d_b, returned_body_data_Size);
	//cudaMalloc((void**)&d_c, double_data_size);

	cudaMalloc((void**)&d_vx, double_data_size);
	cudaMalloc((void**)&d_vy, double_data_size);
	cudaMalloc((void**)&d_rx, double_data_size);
	cudaMalloc((void**)&d_ry, double_data_size);
	cudaMalloc((void**)&d_m, double_data_size);
	cudaMalloc((void**)&d_r_rx, double_data_size);
	cudaMalloc((void**)&d_r_ry, double_data_size);

	//cudaMalloc((void**)&d_bs, body_struct_data_size);


	//vx = (double*)malloc(double_data_size);
	//vy = (double*)malloc(double_data_size);

	//cudaMalloc((void**)&buffer_C, double_data_size);

	int nBlocks = N / BLOCK_SIZE;
	//int nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
	std::cout << "Number of blocks: " << nBlocks << ", block size: " << BLOCK_SIZE << std::endl;
	// Get the start time
	auto start = std::chrono::system_clock::now();
	double dt = 1e11;
	for (int sim_iterations = 0; sim_iterations <= 5000; sim_iterations++)
	{
		//addforces(N);
		//cudaMemcpy(d_a, &bodies[0], body_data_size, cudaMemcpyHostToDevice);

		cudaMemcpy(d_vx, &vx[0], double_data_size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_vy, &vy[0], double_data_size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_rx, &rx[0], double_data_size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_ry, &ry[0], double_data_size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_m, &mass[0], double_data_size, cudaMemcpyHostToDevice);
		//cudaMemcpy(d_bs, &bs[0], body_struct_data_size, cudaMemcpyHostToDevice);
		// Max 1024 threads per block
		
		add_force_simple_doubles<<<nBlocks, BLOCK_SIZE>>>(d_vx, d_vy, d_rx, d_ry, d_m, d_r_rx, d_r_ry, G, N, dt);
		//add_force_simple_struct << <nBlocks, BLOCK_SIZE >> >(d_bs, d_r_rx, d_r_ry, G, N, dt);
		
		//kernel << <1, 1 >> > ();
		// Wait for kernel to complete
		cudaDeviceSynchronize();
		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess)
		{
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(error));
			exit(-1);
		}
		// Read output buffers back to the host
		//cudaMemcpy(f[0], d_b, force_data_size, cudaMemcpyDeviceToHost);
		//cudaMemcpy(&returned_bodies[0], d_b, returned_body_data_Size, cudaMemcpyDeviceToHost);
		//cudaMemcpy(&vx[0], d_b, double_data_size, cudaMemcpyDeviceToHost);
		//cudaMemcpy(&vy[0], d_c, double_data_size, cudaMemcpyDeviceToHost);
		cudaMemcpy(&r_rx[0], d_r_rx, double_data_size, cudaMemcpyDeviceToHost);
		cudaMemcpy(&r_ry[0], d_r_ry, double_data_size, cudaMemcpyDeviceToHost);

		//cudaMemcpy(&fy[0], buffer_C, double_data_size, cudaMemcpyDeviceToHost);
		error = cudaGetLastError();
		if (error != cudaSuccess)
		{
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(error));
			exit(-1);
		}
		for (int i = 0; i < N; i++) 
		{
			//bodies[i]->fx = f[i]->fx;
			//bodies[i]->fy = f[i]->fy;
			//returned_bodies[i].update(1e11);
			//bodies[i]->update(1e11);
			bodies[i]->rx = r_rx[i];
			bodies[i]->ry = r_ry[i];
			// Integrate
			//bodies[i]->rx += dt * vx[i];
			//bodies[i]->ry += dt * vy[i];
		}

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
	auto end = std::chrono::system_clock::now();
	// Get the total time
	auto total = end - start;

	cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(total).count() << " ms" << endl;

	al_destroy_display(display);

	// Clean up resources
	cudaFree(d_vx);
	cudaFree(d_vy);
	cudaFree(d_rx);
	cudaFree(d_ry);
	cudaFree(d_m);
	//cudaFree(d_bs);

	return 0;
}