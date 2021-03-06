#pragma once

#include "cgp/cgp.hpp"

#include "cloth/cloth.hpp"
#include "shape/custom_shape.hpp"
#include "simulation/simulation.hpp"
#include "uniform_grid/Rgrid.hpp"





// The element of the GUI that are not already stored in other structures
struct gui_parameters {
	bool display_frame     = false;
	bool display_wireframe = false;
};



// The structure of the custom scene
struct scene_structure {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	enum DemoScene { sc_CURTAIN, sc_LINEAR, sc_QUADRATIC, sc_TRAMPOLINE };

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
	scene_parameters scene_parameters;
	simulation_parameters parameters;          // Stores the parameters of the simulation (stiffness, mass, damping, time step, etc)
	constraint_structure constraint;           // Handle the parameters of the constraints (fixed vertices, floor and sphere)

	// Particle system related structures
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle
	cgp::buffer<particle_element> all_particles = {};
	cgp::buffer<shape_structure> all_shapes = {};

	// Helper variables
	bool simulation_running = false;   // Boolean indicating if the simulation should be computed
	GLuint cloth_texture;             // Storage of the texture ID used for the cloth


	

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void setShapes(DemoScene demoType = sc_LINEAR);  // Creates the shapes used in the simulation
	void simulate();  // 
	void display(double elapsedTime);     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop
	
	void addCube(float c_x, float c_y, float c_z, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
	
	void addCubeQuadratic(float c_x, float c_y, float c_z, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
	
	void addSphere(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
	
	void addPyramid(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler);
	
	// type 0 = curtain two fixed points
	// type 1 = curtain fixed upper line
	// type 2 = trampoline	
	void addCloth(float c_w, float c_h, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler, int type);
};





