#include "shape.hpp"

using namespace cgp;


void shape_structure::initialize(float c)
{
	// Initial particle spacing (relative to h)
	float const h = p_parameters.h;
	// Fill a square with particles
	particles.clear();
	for (float x = h; x < 0.5f - h; x = x + c * h)
	{
		for (float y = h; y < 0.5f - h; y = y + c * h)
		{
			for (float z = h; z < 0.5f - h; z = z + c * h)
			{
				float x_p = x - 1;
				float y_p = y - 1;
				float z_p = z;
				particle_element particle;
				// rand_interval()
				particle.position = { x_p + h / 8.0, y_p + h / 8.0, z_p + h / 8.0 }; // a zero value in z position will lead to a 2D simulation
				particles.push_back(particle);
			}
		}
	}
}