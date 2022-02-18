#include "simulation.hpp"

using namespace cgp;




// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters)
{
	// direct access to the variables
	grid_2D<vec3>& force = cloth.force;
	grid_2D<vec3> const& position = cloth.position;
	grid_2D<vec3> const& velocity = cloth.velocity;

	size_t const N = cloth.position.size();       // total number of vertices
	size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid

	// Retrive simulation parameter
	float const K = parameters.K;              // spring stifness
	float const m = parameters.mass_total / N; // mass of a particle
	float const mu = parameters.mu;            // damping coefficient
	float const	L0 = 1.0f / (N_edge - 1.0f);   // rest length between two direct neighboring particle


	// Gravity
	const vec3 g = { 0,0,-9.81f };
	for (int ku = 0; ku < N_edge; ++ku)
		for (int kv = 0; kv < N_edge; ++kv)
			force(ku, kv) = m * g;

	// Drag
	for (int ku = 0; ku < N_edge; ++ku)
		for (int kv = 0; kv < N_edge; ++kv)
			force(ku, kv) += -mu * m * velocity(ku, kv);


	// TO DO: Add spring forces ...
	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			for (int du = -2; du <= 2; du++) {
				for (int dv = -2; dv <= 2; dv++) {
					if (std::abs(du) + std::abs(dv) > 2 || std::abs(du) + std::abs(dv) == 0)
						continue;
					int n_ku = ku + du;
					int n_kv = kv + dv;
					if (n_ku >= 0 && n_ku < N_edge && n_kv >= 0 && n_kv < N_edge) {
						vec3 diff = position(n_ku, n_kv) - position(ku, kv);
						float diff_norm = norm(diff);
						vec3 diff_normalized = normalize(diff);
						force(ku, kv) += (K * (diff_norm - (sqrt(std::abs(du) * std::abs(du) + std::abs(dv) * std::abs(dv))) * L0)) * diff_normalized;
					}
				}
			}
		}
	}

	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			force(ku, kv) += parameters.wind.magnitude * dot(cloth.normal(ku, kv), parameters.wind.direction) * m * cloth.normal(ku, kv);
		}
	}
}

void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt)
{
	int const N_edge = cloth.N_samples_edge();
	float const m = parameters.mass_total / static_cast<float>(N_edge);

	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			vec3& v = cloth.velocity(ku, kv);
			vec3& p = cloth.position(ku, kv);
			vec3 const& f = cloth.force(ku, kv);

			// Standard semi-implicit numerical integration
			v = v + dt * f / m;
			p = p + dt * v;
		}
	}

}

void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint)
{
	size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid

	// Fixed positions of the cloth
	for (auto const& it : constraint.fixed_sample) {
		position_contraint c = it.second;
		cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
	}

	// To do: apply external constraints
	// For all vertex:
	//   If vertex is below floor level ...
	//   If vertex is inside collision sphere ...
	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			vec3& v = cloth.velocity(ku, kv);
			vec3& p = cloth.position(ku, kv);
			if (norm(constraint.sphere.center - p) < constraint.sphere.radius + 0.01) {
				vec3 dir = normalize(p - constraint.sphere.center);
				p = constraint.sphere.center + (constraint.sphere.radius + 0.01) * dir;
				vec3 vpar = dot(v, dir) * dir;
				vec3 vn = v - vpar;
				v = vn;
			}
			else if (p.z < constraint.ground_z) {
				p.z = constraint.ground_z + 0.0005;

				vec3 vn = { 0,0,v.z };
				vec3 vplanar = { v.x,v.y,0 };
				v = vplanar;

			}
		}
	}
}

bool simulation_detect_divergence(cloth_structure const& cloth)
{
	bool simulation_diverged = false;
	const size_t N = cloth.position.size();
	for (size_t k = 0; simulation_diverged == false && k < N; ++k)
	{
		const float f = norm(cloth.force.data.at_unsafe(k));
		const vec3& p = cloth.position.data.at_unsafe(k);

		if (std::isnan(f)) // detect NaN in force
		{
			std::cout << "\n **** NaN detected in forces" << std::endl;
			simulation_diverged = true;
		}

		if (f > 600.0f) // detect strong force magnitude
		{
			std::cout << "\n **** Warning : Strong force magnitude detected " << f << " at vertex " << k << " ****" << std::endl;
			simulation_diverged = true;
		}

		if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) // detect NaN in position
		{
			std::cout << "\n **** NaN detected in positions" << std::endl;
			simulation_diverged = true;
		}
	}

	return simulation_diverged;
}

/// <summary>
/// simulation of shape
/// </summary>
/// <param name="shape"></param>
/// <param name="parameters"></param>

void simulation_compute_force(shape_structure& shape, simulation_parameters const& parameters)
{
	size_t const N = shape.particles.size();       // total number of vertices

	// Retrive simulation parameter
	float const K = parameters.K;              // spring stifness
	float const m = parameters.mass_total / N; // mass of a particle
	float const mu = parameters.mu;            // damping coefficient
	//float const	L0 = 1.0f / (N_edge - 1.0f);   // rest length between two direct neighboring particle


	// Gravity
	const vec3 g = { 0,0,-9.81f };
	for (int pIndex = 0; pIndex < N; ++pIndex)
		shape.particles[pIndex].force = m * g;

	/*// Drag
	for (int ku = 0; ku < N_edge; ++ku)
		for (int kv = 0; kv < N_edge; ++kv)
			force(ku, kv) += -mu * m * velocity(ku, kv);


	// TO DO: Add spring forces ...
	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			for (int du = -2; du <= 2; du++) {
				for (int dv = -2; dv <= 2; dv++) {
					if (std::abs(du) + std::abs(dv) > 2 || std::abs(du) + std::abs(dv) == 0)
						continue;
					int n_ku = ku + du;
					int n_kv = kv + dv;
					if (n_ku >= 0 && n_ku < N_edge && n_kv >= 0 && n_kv < N_edge) {
						vec3 diff = position(n_ku, n_kv) - position(ku, kv);
						float diff_norm = norm(diff);
						vec3 diff_normalized = normalize(diff);
						force(ku, kv) += (K * (diff_norm - (sqrt(std::abs(du) * std::abs(du) + std::abs(dv) * std::abs(dv))) * L0)) * diff_normalized;
					}
				}
			}
		}
	}

	for (int ku = 0; ku < N_edge; ++ku) {
		for (int kv = 0; kv < N_edge; ++kv) {
			force(ku, kv) += parameters.wind.magnitude * dot(cloth.normal(ku, kv), parameters.wind.direction) * m * cloth.normal(ku, kv);
		}
	}*/
}

void simulation_numerical_integration(shape_structure& shape, simulation_parameters const& parameters, float dt)
{
	int const N = shape.particles.size();
	float const m = parameters.mass_total / static_cast<float>(N);

	for (int pIndex = 0; pIndex < N; ++pIndex) {

		vec3& v = shape.particles[pIndex].velocity;
		vec3& p = shape.particles[pIndex].position;
		vec3 const& f = shape.particles[pIndex].force;

		// Standard semi-implicit numerical integration
		v = v + dt * f / m;
		p = p + dt * v;
	}

}

void simulation_apply_constraints(shape_structure& shape, constraint_structure const& constraint)
{
	int const N = shape.particles.size();
	/*
	// Fixed positions of the cloth
	for (auto const& it : constraint.fixed_sample) {
		position_contraint c = it.second;
		cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
	}
	*/
	// To do: apply external constraints
	// For all vertex:
	//   If vertex is below floor level ...
	//   If vertex is inside collision sphere ...
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		vec3& v = shape.particles[pIndex].velocity;
		vec3& p = shape.particles[pIndex].position;
		if (norm(constraint.sphere.center - p) < constraint.sphere.radius + 0.01) {
			vec3 dir = normalize(p - constraint.sphere.center);
			p = constraint.sphere.center + (constraint.sphere.radius + 0.01) * dir;
			vec3 vpar = dot(v, dir) * dir;
			vec3 vn = v - vpar;
			v = vn;
		}
		else if (p.z < constraint.ground_z) {
			p.z = constraint.ground_z + 0.0005;

			vec3 vn = { 0,0,v.z };
			vec3 vplanar = { v.x,v.y,0 };
			v = vplanar;

		}
	}

}


/*
void simulate(float dt, buffer<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{

	// Update values
	update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
	update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
	update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
	int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}


	// Collision
	float const epsilon = 1e-3f;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;

		// small perturbation to avoid alignment
		if( p.y<-1 ) {p.y = -1+epsilon*rand_interval();  v.y *= -0.5f;}
		if( p.x<-1 ) {p.x = -1+epsilon*rand_interval();  v.x *= -0.5f;}
		if( p.x>1 )  {p.x =  1-epsilon*rand_interval();  v.x *= -0.5f;}
	}

}
*/