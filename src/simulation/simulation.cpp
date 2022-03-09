#include "simulation.hpp"

using namespace cgp;




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

void simulation_compute_force(cgp::buffer<particle_element>& all_particles, simulation_parameters const& parameters)
{
	size_t const N = all_particles.size();
	// Gravity
	const vec3 g = { 0,0,-9.81f };
	for (int pIndex = 0; pIndex < N; ++pIndex) {
			all_particles[pIndex].force = all_particles[pIndex].mass * g;
	}

}

void calculateOptimalRotation(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	
	
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	

	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		shape.Apq = cgp::mat3();
	}

	for (int i = 0; i < all_particles.size();i++) {
		const particle_element& particle = all_particles[i];
		shape_structure& shape = all_shapes[particle.phase];
		p(0, 0) = particle.position.x-shape.com.x;
		p(1, 0) = particle.position.y-shape.com.y;
		p(2, 0) = particle.position.z-shape.com.z;

		q(0, 0) = shape.relativeLocations[i-shape.relativeLocationsOffset].x;
		q(0, 1) = shape.relativeLocations[i-shape.relativeLocationsOffset].y;
		q(0, 2) = shape.relativeLocations[i-shape.relativeLocationsOffset].z;

		
		
		shape.Apq += particle.mass * p * q;
	}


	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		
		cgp::mat3 sym = transpose(shape.Apq) * shape.Apq;
		Eigen::MatrixXd A(3, 3);
		A(0, 0) = sym(0, 0);
		A(0, 1) = sym(0, 1);
		A(0, 2) = sym(0, 2);

		A(1, 0) = sym(1, 0);
		A(1, 1) = sym(1, 1);
		A(1, 2) = sym(1, 2);

		A(2, 0) = sym(2, 0);
		A(2, 1) = sym(2, 1);
		A(2, 2) = sym(2, 2);


		Eigen::EigenSolver <Eigen::MatrixXd > es(A);
		Eigen::MatrixXcd Dc = es.eigenvalues().asDiagonal();
		Eigen::MatrixXcd Vc = es.eigenvectors();

		Eigen::MatrixXd D = Dc.real();
		Eigen::MatrixXd V = Vc.real();

		D(0, 0) = 1.f / sqrt(D(0, 0));

		D(1, 1) = 1.f / sqrt(D(1, 1));
		D(2, 2) = 1.f / sqrt(D(2, 2));

		Eigen::MatrixXd Sinv = (V * D * V.inverse());

		cgp::mat3 S_ = cgp::mat3();

		S_(0, 0) = Sinv(0, 0);
		S_(0, 1) = Sinv(0, 1);
		S_(0, 2) = Sinv(0, 2);


		S_(1, 0) = Sinv(1, 0);
		S_(1, 1) = Sinv(1, 1);
		S_(1, 2) = Sinv(1, 2);


		S_(2, 0) = Sinv(2, 0);
		S_(2, 1) = Sinv(2, 1);
		S_(2, 2) = Sinv(2, 2);

		cgp::mat3 result = shape.Apq * S_;

		if (det(result) < 0) {
			result = result / det(result);
		}
		//TODO print determinant

		shape.optimalRotation = result;
		shape.A = shape.Apq * shape.Aqq;

		shape.A = shape.A / pow(cgp::det(shape.A), 1.f / 3.f);
		//std::cout << shape.optimalRotation << std::endl;
	}
	

	

}

void calculateCurrentCom(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	
	for (int i = 0; i < all_shapes.size(); i++) {
		all_shapes[i].com = vec3(0.0, 0.0, 0.0);
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		all_shapes[particle.phase].com += particle.position;
	}

	for (int i = 0; i < all_shapes.size(); i++) {
		all_shapes[i].com /= all_shapes[i].nbParticles;
	}


}

void preCalculations(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	calculateCurrentCom(all_particles,all_shapes);
	calculateOptimalRotation(all_particles, all_shapes);
	
}

void simulation_numerical_integration(cgp::buffer<particle_element>& all_particles, float dt)
{

	int const N = all_particles.size();
	for (int pIndex = 0; pIndex < N; ++pIndex) {

		vec3& v = all_particles[pIndex].velocity;
		vec3& x = all_particles[pIndex].position;


		
		vec3 const& f = all_particles[pIndex].force;
		float m = all_particles[pIndex].mass;
		v = v + dt * f / m ; //simplification needed
		x = x + dt * v;
	}

}

void shapeMatching(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	const float alpha = 0.8;
	int const N = all_particles.size();
	
	cgp::buffer<vec3> g = {}; //size N at the end
	
	cgp::buffer<cgp::vec3> coms;
	coms.resize(all_shapes.size());
	coms.fill(vec3(0.0, 0.0, 0.0));

	for (int pIndex = 0; pIndex < N; ++pIndex) {
		const particle_element& particle = all_particles[pIndex];
		shape_structure& shape = all_shapes[particle.phase];
		vec3 res = (shape.optimalRotation * 0.2f + shape.A * 0.8f) * (shape.relativeLocations[pIndex-shape.relativeLocationsOffset]);
		
		g.push_back(res);
		coms[particle.phase] += res;
	}
	
	for (int i = 0; i < all_shapes.size(); i++) {
		coms[i] /= all_shapes[i].nbParticles;
		std::cout << coms[i] << std::endl;
	}

	for (int pIndex = 0; pIndex < N; ++pIndex) {
		const particle_element& particle = all_particles[pIndex];
		g(pIndex) = g(pIndex) - coms[particle.phase];
	}
	
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		particle_element& particle = all_particles[pIndex];
		shape_structure& shape = all_shapes[particle.phase];
		vec3& v = particle.velocity;
		vec3& x = particle.position;
		vec3 target = g[pIndex] + shape.com;
		x = x + alpha * (target - x);
	}
}


void simulation_apply_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float di)
{
	int const N = all_particles.size();
	const vec3 g = { 0,0,-9.81f };
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		vec3& v = all_particles[pIndex].velocity;
		vec3& p = all_particles[pIndex].position;
		for (const auto& wall : constraint.walls) {
			
			if (dot(p - wall.point, wall.normal) < 0.01) {
				vec3 dp = (-dot(p - wall.point, wall.normal) + 0.01) * wall.normal;
				p += dp;
				prevX[pIndex] += dp;
				all_particles[pIndex].flagConstraint = true;
				vec3 vn = dot(v, wall.normal) * wall.normal;
				vec3 vpar = v - vn;
				all_particles[pIndex].dv = (0.5* vpar - vn) - v;
			}
		}
		for (const auto& sphere : constraint.spheres) {

			if (norm(sphere.center - p) < sphere.radius + 0.01) {
				vec3 dir = normalize(p - sphere.center);
				p = sphere.center + (sphere.radius + 0.01) * dir;
				all_particles[pIndex].flagConstraint = true;
				vec3 vpar = dot(v, dir) * dir;
				vec3 vn = v - vpar;
				all_particles[pIndex].dv = (vn - vpar) - v;
			}
		}
		
	

	}
		
}


void adjustVelocity(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, float dt){
	for (int index = 0; index < all_particles.size(); index++) {
		particle_element& particle = all_particles[index];
		particle.velocity = (particle.position - prevX[index]) / dt;
		if (particle.flagConstraint) {
			particle.flagConstraint = false;
			particle.velocity += particle.dv;
			particle.dv = vec3(0,0,0);
		}
		
		if (norm(particle.position - prevX[index]) < 0.001) {
			particle.position = prevX[index];
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