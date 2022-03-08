#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"

// Parameters of the colliding sphere (center, radius)
struct sphere_parameter {
	cgp::vec3 center;
	float radius;
};

// Parameter attached to a fixed vertex (ku,kv) coordinates + 3D position
struct position_contraint {
	int ku;
	int kv;
	cgp::vec3 position;
};

struct plane_contraint {
	cgp::vec3 point;
	cgp::vec3 normal;
};


//cloth constraint_structure
struct constraint_structure
{
	float cubeSize = 1.5;
	cgp::buffer<plane_contraint> walls = { { {0.0f, 0.0f, 0.0f},  {0.0f, 0.0f, 1.0f}}, 
											{ {cubeSize, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}},
											{ {-cubeSize, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
											{ {0.0f, cubeSize, 0.0f}, {0.0f, -1.0f, 0.0f}},
											{ {0.0f, -cubeSize, 0.0f}, {0.0f, 1.0f, 0.0f}} };                                   // Height of the flood
	cgp::buffer<sphere_parameter> spheres = { {{0.1f, 0.5f, 0.2f}, 0.15f},
											{ {-0.7f, 0.5f, 0.2f}, 0.15f},
											{ {-0.7f, 1.25f, 0.4f}, 0.15f},
	}; // Colliding sphere
	
	std::map<size_t, position_contraint> fixed_sample; // Storage of all fixed position of the cloth

	// Add a new fixed position
	void add_fixed_position(int ku, int kv, cloth_structure const& cloth);
	// Remove a fixed position
	void remove_fixed_position(int ku, int kv);

};
