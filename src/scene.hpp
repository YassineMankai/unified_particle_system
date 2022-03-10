#pragma once

#include "cgp/cgp.hpp"

#include "cloth/cloth.hpp"
#include "shape/shape.hpp"
#include "simulation/simulation.hpp"

// The element of the GUI that are not already stored in other structures
struct gui_parameters {
	bool display_frame     = false;
	bool display_wireframe = false;
	bool display_particles = true;
	bool display_color = true;
	int N_sample_edge    = 20;  // number of samples of the cloth (the total number of vertices is N_sample_edge^2)
};



// The structure of the custom scene
struct scene_structure {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	cgp::mesh_drawable global_frame;          // The standard global frame
	cgp::scene_environment_basic environment; // Standard environment controler
	gui_parameters gui;                       // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	// Display of the obstacles and constraints
	cgp::mesh_drawable obstacle_floor;
	cgp::mesh_drawable obstacle_sphere;

	// simulation related structures
	particle_parameters_structure p_parameters;
	simulation_parameters parameters;          // Stores the parameters of the simulation (stiffness, mass, damping, time step, etc)
	constraint_structure constraint;           // Handle the parameters of the constraints (fixed vertices, floor and sphere)

	// Particle system related structures
	shape_structure shape;
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle

	/*
	cgp::grid_2D<cgp::vec3> field;      // grid used to represent the volume of the fluid under the particles
	cgp::mesh_drawable field_quad; // quad used to display this field color
	*/

	// Helper variables
	bool simulation_running = true;   // Boolean indicating if the simulation should be computed
	GLuint cloth_texture;             // Storage of the texture ID used for the cloth


	cgp::buffer<particle_element> all_particles = {};
	cgp::buffer<shape_structure> all_shapes = {};

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void display(double elapsedTime);     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop
	void addCube(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
	void addCubeQuadratic(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
};





