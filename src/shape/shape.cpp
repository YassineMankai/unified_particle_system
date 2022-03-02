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
				particle.position =vec3( x_p + h / 8.0f, y_p + h / 8.0f, z_p + h / 8.0f)+vec3(0.0,0.0,1.0f); // a zero value in z position will lead to a 2D simulation
				centerOfMass += particle.position;
				sum++;
				particles.push_back(particle);
			}
		}
	}
	//initializing center of mass
	centerOfMass /= sum;
	com0 = centerOfMass;
	
	
	rotation_transform rotation = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), 75);
	affine_rt rt = rotation_around_center(rotation, com0);


	for (int i = 0; i < particles.size(); i++) {
		particles[i].position = (rotation * (particles[i].position - com0)) + com0;
	};

	vec3 res = vec3(0.0, 0.0, 0.0);
	for (const particle_element& particle : particles) {
		res = res + particle.position;
	}
	res = res / sum;

	if (norm(res - com0) < 0.001) {
		for (particle_element& particle : particles) {
			particle.position = particle.position - res + com0;
		}
	}
	

	//initializing relative locations
	relativeLocations = {};
	for (int i = 0; i < particles.size();i++) {
		const particle_element& particleP = particles[i];
		relativeLocations.push_back(particleP.position - com0);
	}
	
	//pre calculation of Aqq
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	for (int i = 0; i < particles.size(); i++) {
		const particle_element& particle = particles[i];
		p(0, 0) = relativeLocations[i].x;
		p(1, 0) = relativeLocations[i].y;
		p(2, 0) = relativeLocations[i].z;

		q(0, 0) = relativeLocations[i].x;
		q(0, 1) = relativeLocations[i].y;
		q(0, 2) = relativeLocations[i].z;

		Aqq += particle.mass * p * q;
	}



}