#include "simulation.hpp"

using namespace cgp;



int phaseToShapeIndex(int phase) {
	return std::abs(phase) - 1;
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

void simulation_compute_force(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, simulation_parameters const& parameters)
{
	size_t const N = all_particles.size();
	// Gravity
	const vec3 g = { 0,0,-9.8f };
	for (int pIndex = 0; pIndex < N; ++pIndex) {
			all_particles[pIndex].force = all_particles[pIndex].mass * g;
		
	}
}

void calculateOptimalRotation(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();

	cgp::mat39 pQuad = cgp::mat39();
	cgp::mat9 qQuad = cgp::mat9();


	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		if (shape.isShapeMatching) {
			if (!shape.isQuadratic)
				shape.Apq = cgp::mat3();
			else
			{
				shape.Apq = cgp::mat3();
				shape.ApqQuad = cgp::mat39();
			}
		}
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		if (particle.phase >= 0) {
			shape_structure& shape = all_shapes[phaseToShapeIndex(particle.phase)];
			if (!shape.isQuadratic) {
				p(0, 0) = particle.position.x - shape.com.x;
				p(1, 0) = particle.position.y - shape.com.y;
				p(2, 0) = particle.position.z - shape.com.z;

				q(0, 0) = shape.relativeLocations[i - shape.relativeLocationsOffset].x;
				q(0, 1) = shape.relativeLocations[i - shape.relativeLocationsOffset].y;
				q(0, 2) = shape.relativeLocations[i - shape.relativeLocationsOffset].z;

				shape.Apq += particle.mass * p * q;
			}
			else {
				p(0, 0) = particle.position.x - shape.com.x;
				p(1, 0) = particle.position.y - shape.com.y;
				p(2, 0) = particle.position.z - shape.com.z;

				q(0, 0) = shape.relativeLocations[i - shape.relativeLocationsOffset].x;
				q(0, 1) = shape.relativeLocations[i - shape.relativeLocationsOffset].y;
				q(0, 2) = shape.relativeLocations[i - shape.relativeLocationsOffset].z;

				shape.Apq += particle.mass * p * q;

				cgp::vec9& q = shape.qQuad[i - shape.relativeLocationsOffset];

				pQuad(0, 0) = particle.position.x - shape.com.x;
				pQuad(1, 0) = particle.position.y - shape.com.y;
				pQuad(2, 0) = particle.position.z - shape.com.z;

				qQuad(0, 0) = q(0);
				qQuad(0, 1) = q(1);
				qQuad(0, 2) = q(2);

				qQuad(0, 3) = q(3);
				qQuad(0, 4) = q(4);
				qQuad(0, 5) = q(5);

				qQuad(0, 6) = q(6);
				qQuad(0, 7) = q(7);
				qQuad(0, 8) = q(8);


				shape.ApqQuad += particle.mass * pQuad * qQuad;
			}
		}
	}
	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		if (shape.isShapeMatching) {
			if (!shape.isQuadratic) {
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

				shape.optimalRotation = result;
				shape.A = shape.Apq * shape.Aqq;
				shape.A = shape.A / pow(cgp::det(shape.A), 1.f / 3.f);
			}
			else {
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

				cgp::mat39& rotation = shape.optimalRotationQuad;

				rotation = cgp::mat39();

				rotation(0, 0) = result(0, 0);
				rotation(0, 1) = result(0, 1);
				rotation(0, 2) = result(0, 2);

				rotation(1, 0) = result(1, 0);
				rotation(1, 1) = result(1, 1);
				rotation(1, 2) = result(1, 2);

				rotation(2, 0) = result(2, 0);
				rotation(2, 1) = result(2, 1);
				rotation(2, 2) = result(2, 2);


				cgp::mat39 aQuad = shape.ApqQuad * shape.AqqQuad;

				cgp::mat3 A_ = cgp::mat3();

				A_(0, 0) = aQuad(0, 0);
				A_(1, 0) = aQuad(1, 0);
				A_(2, 0) = aQuad(2, 0);

				A_(0, 1) = aQuad(0, 1);
				A_(1, 1) = aQuad(1, 1);
				A_(2, 1) = aQuad(2, 1);

				A_(0, 2) = aQuad(0, 2);
				A_(1, 2) = aQuad(1, 2);
				A_(2, 2) = aQuad(2, 2);

				A_ = A_ / pow(cgp::det(A_), 1.f / 3.f);

				aQuad(0, 0) = A_(0, 0);
				aQuad(1, 0) = A_(1, 0);
				aQuad(2, 0) = A_(2, 0);

				aQuad(0, 1) = A_(0, 1);
				aQuad(1, 1) = A_(1, 1);
				aQuad(2, 1) = A_(2, 1);

				aQuad(0, 2) = A_(0, 2);
				aQuad(1, 2) = A_(1, 2);
				aQuad(2, 2) = A_(2, 2);

				shape.AQuad = aQuad;



			}
		}

	}
}

void calculateCurrentCom(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	
	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].isShapeMatching)
		all_shapes[i].com = vec3(0.0, 0.0, 0.0);
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		if (particle.phase>=0)
		all_shapes[phaseToShapeIndex(particle.phase)].com += particle.position;
	}

	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].isShapeMatching)
		all_shapes[i].com /= all_shapes[i].nbParticles;
	}


}

void preCalculations(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	calculateCurrentCom(all_particles, all_shapes);
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
		v = v + dt * f / m; //simplification needed
		x = x + dt * v;
	}

}

void shapeMatching(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, float alpha, float beta) {

	//const float alpha = 0.3;
	int const N = all_particles.size();

	cgp::buffer<vec3> g = {}; //size N at the end
	
	cgp::buffer<cgp::vec3> coms;
	int numberOfShapeMatching = 0;

	coms.resize(all_shapes.size());
	coms.fill(cgp::vec3(0.0, 0.0, 0.0));

	for (int pIndex = 0; pIndex < N; ++pIndex) {
		const particle_element& particle = all_particles[pIndex];
		vec3 res = vec3(0.0, 0.0, 0.0);
		if (particle.phase >= 0) {
			shape_structure& shape = all_shapes[phaseToShapeIndex(particle.phase)];
			
			if (!shape.isQuadratic) {
				res = (shape.optimalRotation * beta + shape.A * (1 - beta)) * (shape.relativeLocations[pIndex - shape.relativeLocationsOffset]);
			}
			else {
				res = ((1 - beta) * shape.AQuad + beta * shape.optimalRotationQuad) * shape.qQuad[pIndex - shape.relativeLocationsOffset];
			}
			
			coms[phaseToShapeIndex(particle.phase)] += res;
		}
		g.push_back(res); 
		//can be optimised to be filled only with targets for shape matching for now if the particle doesn't belong to 
		//shape matching object the res stays 0
	}
	
	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].isShapeMatching)
		coms[i] /= all_shapes[i].nbParticles;
	}
	
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		const particle_element& particle = all_particles[pIndex];
		if (particle.phase >= 0) {
			g(pIndex) = g(pIndex) - coms[phaseToShapeIndex(particle.phase)];
		}
	}
	
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		particle_element& particle = all_particles[pIndex];
		if (particle.phase >= 0) {
			shape_structure& shape = all_shapes[phaseToShapeIndex(particle.phase)];
			vec3& x = particle.position;
			vec3 target = g[pIndex] + shape.com;
			x = x + alpha * (target - x);
		}
	}
	
}

void simulation_apply_shape_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, constraint_structure const& constraint) {
	 


	
	for (int pIndex = 0; pIndex < all_particles.size(); pIndex++) {
		particle_element& particle = all_particles[pIndex];
		if (particle.phase < 0) { //not shapematching
			
			shape_structure& shape = all_shapes[phaseToShapeIndex(particle.phase)];
			vec3& x = particle.position;
			float k = 3.0; //stiffness parameter TODO put it in the parameters of simulation
			int index1, index2, index3, index4;
			float L;
			//int numberOfConstraints = shape.width*shape.height*3; //204
			int numberOfConstraints = 12; //204
			//structure constraints
			L = shape.structLength0;
			index1 = shape.getNeighbour(pIndex, 1, 0);		
			index2 = shape.getNeighbour(pIndex, -1, 0);
			index3 = shape.getNeighbour(pIndex, 0, 1);
			index4 = shape.getNeighbour(pIndex, 0, -1);


			if (index1 != -1) {
				const particle_element& p = all_particles[index1];
				float  distance = norm(p.position - particle.position);
				x += k *(distance - L) * normalize(p.position - particle.position)/numberOfConstraints;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				x += k  * (distance - L) * normalize(p.position - particle.position)/numberOfConstraints;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position)/numberOfConstraints;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				x += k  * (distance - L) * normalize(p.position - particle.position)/numberOfConstraints;
			}
			//shear constraints
			L = shape.shearLength0;
			index1 = shape.getNeighbour(pIndex, 1, 1);
			index2 = shape.getNeighbour(pIndex, -1, -1);
			index3 = shape.getNeighbour(pIndex, 1, -1);
			index4 = shape.getNeighbour(pIndex, -1, 1);


			if (index1 != -1) {
				const particle_element& p = all_particles[index1];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}


			//bending constraints
			L = shape.bendLength0;
			index1 = shape.getNeighbour(pIndex, 0, 2);

			index2 = shape.getNeighbour(pIndex, 0, -2);
			index3 = shape.getNeighbour(pIndex, 2, 0);
			index4 = shape.getNeighbour(pIndex, -2, 0);


			if (index1 != -1) {
				const particle_element& p = all_particles[index1];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				x += k * (distance - L) * normalize(p.position - particle.position) / numberOfConstraints;
			}

		}
	}
	
	//fixed points of the scene:
	for (auto const& it : constraint.fixed_sample) {
		int pIndex = it.first;
		particle_element& particle = all_particles[pIndex];
		particle.position = it.second;
	}

	

}


void simulation_apply_env_contact_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float di)
{
	int const N = all_particles.size();

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
				//all_particles[pIndex].dv = (0.5* vpar - 0.7*vn) - v; //cube
				all_particles[pIndex].dv = (0.5 * vpar - 0.8 * vn) - v; //sphere
			}
		}
		for (const auto& sphere : constraint.spheres) {

			if (norm(sphere.center - p) < sphere.radius + 0.01) {
				vec3 dir = normalize(p - sphere.center);
				p = sphere.center + (sphere.radius + 0.01) * dir;
				all_particles[pIndex].flagConstraint = true;
				vec3 vpar = dot(v, dir) * dir;
				vec3 vn = v - vpar;
				all_particles[pIndex].dv = (vn)-v;
			}
		}
	}
}

void simulation_apply_particle_contact_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, Rgrid const& grid, float di)
{
	for (int p1 = 0; p1 < all_particles.size(); p1++)
	{
		particle_element& particle1 = all_particles[p1];
		float const epsilon = 1e-5f;
		cgp::buffer<int> neighborhood = grid.getNeighborhood(particle1.position);

		for (int i = 0; i < neighborhood.size(); i++)
		{
			int p2 = neighborhood[i];
			particle_element& particle2 = all_particles[p2];
			if (particle1.phase == particle2.phase)
				continue;
			float detection = norm(particle1.position - particle2.position);
			if (detection <= 0.01f * 2) { //TODO: use correct sphere radius
				vec3 u = (particle1.position - particle2.position) / detection;
				float diff = 0.01f * 2 - detection;
				vec3 dx = (epsilon + diff / 2.0) * u;
				particle1.position = particle1.position + dx;
				prevX[i] += dx;
				particle2.position = particle2.position - dx;
				prevX[p2] -= dx;
				float projected_relative_speed = dot(particle2.velocity - particle1.velocity, u);
				particle1.flagConstraint = true;
				particle2.flagConstraint = true;
				if (std::abs(projected_relative_speed) > 0.1) {
					particle1.dv += projected_relative_speed * u;
					particle2.dv += -projected_relative_speed * u;
				}
				else
				{
					particle1.dv += -particle1.velocity / 3.f;
					particle2.dv += -particle2.velocity / 3.f;
				}

			}
		}
	}
}

void adjustVelocity(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, float dt) {
	for (int index = 0; index < all_particles.size(); index++) {
		particle_element& particle = all_particles[index];
		particle.velocity = (particle.position - prevX[index]) / dt;
		if (particle.flagConstraint) {
			particle.flagConstraint = false;
			particle.velocity += particle.dv;
			particle.dv = vec3(0, 0, 0);
		}

		if (norm(particle.position - prevX[index]) < 0.0005) {
			particle.position = prevX[index];
		}
	}

}


