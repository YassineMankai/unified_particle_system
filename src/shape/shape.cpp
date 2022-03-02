#include "shape.hpp"

using namespace cgp;


void shape_structure::initialize(float c)
{
	// Initial particle spacing (relative to h)
	float const h = p_parameters.h;
	// Fill a square with particles
	particles = {};
	
	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0,0.0,0.0);
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
				//particle.mass = x/10.f;
				// rand_interval()
				particle.position =vec3( x_p + h / 8.0, y_p + h / 8.0, z_p + h / 8.0)+vec3(0.0,0.0,0.5); // a zero value in z position will lead to a 2D simulation
				centerOfMass += particle.position;
				sum++;
				particles.push_back(particle);
			}
		}
	}
	//initializing center of mass
	centerOfMass /= sum;
	com0 = centerOfMass;

	
	rotation_transform rotation = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), 45);
	affine_rt rt = rotation_around_center(rotation, com0);

	for (int i = 0; i < particles.size(); i++) {
		particles[i].position =(rt.matrix()* vec4(particles[i].position,1.0)).xyz();
	}
	






	//initializing relative locations
	relativeLocations = {};
	for (int i = 0; i < particles.size();i++) {
		const particle_element& particleP = particles[i];
		relativeLocations.push_back(particleP.position - com0);
	}

}