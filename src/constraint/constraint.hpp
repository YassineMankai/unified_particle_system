#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"

// Parameters of the colliding sphere (center, radius)
struct sphere_parameter {
	cgp::vec3 center;
	float radius;
};




struct plane_contraint {
	cgp::vec3 point;
	cgp::vec3 normal;
};


//constraint_structure
struct constraint_structure
{
	float cubeSize = 1.5;
	cgp::buffer<plane_contraint> walls = { { {0.0f, 0.0f, 0.0f},  {0.0f, 0.0f, 1.0f}},
											{ {cubeSize, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}},
											{ {-cubeSize, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
											{ {0.0f, cubeSize, 0.0f}, {0.0f, -1.0f, 0.0f}},
											{ {0.0f, -cubeSize, 0.0f}, {0.0f, 1.0f, 0.0f}} };// Height of the flood
	cgp::buffer<sphere_parameter> spheres = {   {{-0.22f, 0.0f, 0.15f}, 0.15f},
												{ {0.22f, 0.0f, 0.15f}, 0.15f},
												{ {0.0f, -0.22f, 0.15f}, 0.15f},
												{ {0.0f, 0.22f, 0.15f}, 0.15f},
	}; // Colliding sphere
	
	std::map<int, cgp::vec3> fixed_sample; // Storage of all fixed position of the cloth
};
