#include "simulation.hpp"

using namespace cgp;



mat9 calculateInverseWithEigen(mat9 A) {
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 9);

	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 9; j++) {
			B(i, j) = A(i, j);
		}
	}

	B = B.inverse();

	mat9 K = cgp::mat9();

	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 9; j++) {
			K(i, j) = B(i, j);
		}
	}

	return K;


}



// ################# Constraints #################

// Helper functions
void calculateOptimalRotation(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();

	cgp::mat39 pQuad = cgp::mat39();
	cgp::mat9 qQuad = cgp::mat9();

	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		shape.Apq = cgp::mat3();
		shape.ApqQuad = cgp::mat39();
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		shape_structure& shape = all_shapes[particle.phase];
		if (shape.type == RIGID) {
			p(0, 0) = particle.position.x - shape.com.x;
			p(1, 0) = particle.position.y - shape.com.y;
			p(2, 0) = particle.position.z - shape.com.z;

			q(0, 0) = shape.relativeLocations[i - shape.relativeLocationsOffset].x;
			q(0, 1) = shape.relativeLocations[i - shape.relativeLocationsOffset].y;
			q(0, 2) = shape.relativeLocations[i - shape.relativeLocationsOffset].z;

			shape.Apq += particle.mass * p * q;
		}
		else if (shape.type == QUADRATIC) {
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

	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];

		if (shape.type == RIGID) {
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
		else if (shape.type == QUADRATIC) {
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
void calculateCurrentCom(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {

	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].type == RIGID || all_shapes[i].type == QUADRATIC)
			all_shapes[i].com = vec3(0.0, 0.0, 0.0);
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		if (all_shapes[particle.phase].type == RIGID || all_shapes[particle.phase].type == QUADRATIC)
			all_shapes[particle.phase].com += particle.position;
	}

	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].type == RIGID || all_shapes[i].type == QUADRATIC)
			all_shapes[i].com /= all_shapes[i].nbParticles;
	}


}
void preCalculations(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes) {
	calculateCurrentCom(all_particles, all_shapes);
	calculateOptimalRotation(all_particles, all_shapes);
}

// Constraints related to shape type (shape matching, cloth, fluid)
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

		shape_structure& shape = all_shapes[particle.phase];

		if (shape.type == RIGID) {
			res = (shape.optimalRotation * beta + shape.A * (1 - beta)) * (shape.relativeLocations[pIndex - shape.relativeLocationsOffset]);
		}
		else if (shape.type == QUADRATIC) {
			res = ((1 - beta) * shape.AQuad + beta * shape.optimalRotationQuad) * shape.qQuad[pIndex - shape.relativeLocationsOffset];
		}

		coms[particle.phase] += res;

		g.push_back(res);
		//can be optimised to be filled only with targets for shape matching for now if the particle doesn't belong to 
		//shape matching object the res stays 0
	}

	for (int i = 0; i < all_shapes.size(); i++) {
		if (all_shapes[i].type == RIGID || all_shapes[i].type == QUADRATIC)
			coms[i] /= all_shapes[i].nbParticles;
	}

	for (int pIndex = 0; pIndex < N; ++pIndex) {
		const particle_element& particle = all_particles[pIndex];
		if (particle.phase >= 0) {
			g(pIndex) = g(pIndex) - coms[particle.phase];
		}
	}

	for (int pIndex = 0; pIndex < N; ++pIndex) {
		particle_element& particle = all_particles[pIndex];
		shape_structure& shape = all_shapes[particle.phase];
		if (shape.type == RIGID || shape.type == QUADRATIC) {
			vec3& x = particle.position;
			vec3 target = g[pIndex] + shape.com;
			x = x + alpha * (target - x);
		}
	}
}
void simulation_apply_shape_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, constraint_structure const& constraint, simulation_parameters const& parameters) {
	for (int pIndex = 0; pIndex < all_particles.size(); pIndex++) {
		particle_element& particle = all_particles[pIndex];
		shape_structure& shape = all_shapes[particle.phase];
		if (shape.type == CLOTH) {
			vec3& x = particle.position;
			vec3 dx = vec3(0.0, 0.0, 0.0);

			float k = parameters.clothStiffness; //stiffness parameter TODO put it in the parameters of simulation
			int index1, index2, index3, index4;
			float L;
			//int numberOfConstraints = shape.width*shape.height*3; //204
			int numberOfConstraints = 0; //204

			//structure constraints
			L = shape.structLength0;
			index1 = shape.getNeighbour(pIndex, 1, 0);
			index2 = shape.getNeighbour(pIndex, -1, 0);
			index3 = shape.getNeighbour(pIndex, 0, 1);
			index4 = shape.getNeighbour(pIndex, 0, -1);


			if (index1 != -1) {
				const particle_element& p = all_particles[index1];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
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
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
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
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index2 != -1) {
				const particle_element& p = all_particles[index2];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index3 != -1) {
				const particle_element& p = all_particles[index3];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			if (index4 != -1) {
				const particle_element& p = all_particles[index4];
				float  distance = norm(p.position - particle.position);
				dx += k * (distance - L) * normalize(p.position - particle.position);
				numberOfConstraints++;
			}
			x += dx / numberOfConstraints;

		}
	}
}

// Constraints related to the environment (walls, obstacles && fixed points)
void simulation_apply_contact_constraints(cgp::buffer<particle_element>& all_particles, const cgp::buffer<shape_structure>& all_shapes, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float dt)
{
	std::vector<Rgrid> regularGrids;
	for (int i = 0; i < all_shapes.size(); i++) {
		regularGrids.push_back(Rgrid());
		regularGrids.back().initialize(20);
	}

	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		regularGrids[particle.phase].updateMinMax(particle.position);
	}
	
	for (int i = 0; i < all_particles.size(); i++) {
		const particle_element& particle = all_particles[i];
		regularGrids[particle.phase].insert(particle.position, i);
	}
	
	// env contacts
	int const N = all_particles.size();
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		vec3& v = all_particles[pIndex].velocity;
		vec3& p = all_particles[pIndex].position;
		vec3& dx = all_particles[pIndex].dx;
		for (const auto& wall : constraint.walls) {

			if (dot(p - wall.point, wall.normal) < 0.01) {
				vec3 dp = (-dot(p - wall.point, wall.normal) + 0.01) * wall.normal;
				dx += dp;
				all_particles[pIndex].nbConstraint += 1;
				vec3 vn = dot(v, wall.normal) * wall.normal;
				vec3 vpar = v - vn;
				all_particles[pIndex].dx_friction_and_restitution += (-0.7 * vn + 0.4* vpar - v) * dt; //cube
				all_particles[pIndex].nbConstraintTest += 1;
			}
		}
		for (const auto& sphere : constraint.spheres) {

			if (norm(sphere.center - p) < sphere.radius + 0.01) {
				vec3 dir = (p - sphere.center) / norm(p - sphere.center);
				dx += sphere.center + (sphere.radius + 0.01) * dir - p;
				all_particles[pIndex].nbConstraint += 1;
				vec3 vn = dot(v, dir) * dir;
				vec3 vpar = v - vn;
				all_particles[pIndex].dx_friction_and_restitution += (-0.7 * vn + 0.4 * vpar - v) * dt; //cube
				all_particles[pIndex].nbConstraintTest += 1;
			}
		}
	}


	// particle contact
	for (int p1 = 0; p1 < all_particles.size(); p1++)
	{
		particle_element& particle1 = all_particles[p1];
		float const epsilon = 1e-5f;
		cgp::buffer<int> neighborhood;
		for (int i = 0; i < all_shapes.size(); i++) {
			neighborhood.push_back(regularGrids[i].getNeighborhood(particle1.position));
		}
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
				particle1.dx += dx;
				particle2.dx += - dx;

				particle1.nbConstraint += 1;
				particle2.nbConstraint += 1;
				particle1.nbConstraintTest += 1;
				particle2.nbConstraintTest += 1;
				
				if (norm(dx) > 0.005 * dt) {
					float projected_correction = 2.6 * dot(dx, u);
					particle1.dx_friction_and_restitution += projected_correction * u;
					particle2.dx_friction_and_restitution += -projected_correction * u;
				}
				else
				{
					particle1.dx_friction_and_restitution += -particle1.velocity * dt / 3.f;
					particle2.dx_friction_and_restitution += -particle2.velocity * dt / 3.f;
				}
				

				//friction
				/*vec3 dx_fr = (particle1.position - prevX[p1]) - (particle2.position - prevX[p2]);
				u = normalize(particle1.position - particle2.position);
				vec3 dx_fr_normal = dot(dx_fr, u) * u;
				vec3 dx_fr_tangential = dx_fr - dx_fr_normal;
				vec3 dx_fr_f;
				if (norm(dx_fr_normal) < 0.57 * 0.01) {
					dx_fr_f = 0.5f * dx_fr_normal;
				}
				else
				{
					dx_fr_f = 0.5f * dx_fr_normal * std::min(1.0f, 0.4f * 0.02f / norm(dx_fr_normal));
				}
				particle1.dx += dx_fr_f;
				particle2.dx -= 0.5f * dx_fr_f;*/
			}
		}
	}

	// Apply dx
	for (int i = 0; i < all_particles.size(); i++) {
		particle_element& particle = all_particles[i];
		if (particle.nbConstraint > 0) {
			particle.position += particle.dx / particle.nbConstraint;
			prevX[i] += particle.dx / particle.nbConstraint;
			if (particle.nbConstraintTest > 0) {
				prevX[i] -= particle.dx_friction_and_restitution / particle.nbConstraintTest;
				particle.nbConstraintTest = 0;
				particle.dx_friction_and_restitution = vec3(0,0,0);
			}
			particle.nbConstraint = 0;
			particle.dx = vec3(0, 0, 0);
		}
	}

}


// ################# Position based dynamics integration #################
void simulation_compute_force(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, simulation_parameters const& parameters)
{
	size_t const N = all_particles.size();
	// Gravity
	const vec3 g = { 0,0,-9.8f };
	int h = parameters.forceField.dimension(0);
	int w = parameters.forceField.dimension(1);
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		all_particles[pIndex].force = all_particles[pIndex].mass * g;
		vec2 pos2DNormalized = (all_particles[pIndex].position.xy() - vec2(-1.6, -1.6)) / 3.4; // TODO : use a more generic formula
		int2 cell = int2(pos2DNormalized.x * (h-1), pos2DNormalized.y * (w-1));
		all_particles[pIndex].force += parameters.forceField[cell];
	}
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
void adjustVelocity(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, float dt) {
	for (int index = 0; index < all_particles.size(); index++) {
		particle_element& particle = all_particles[index];
		particle.velocity = (particle.position - prevX[index]) / dt;
		
		//can consider freezing time per shape, it depends on the simulation
		if (norm(particle.position - prevX[index]) < 0.0001) {
			particle.position = prevX[index];
		}
	}

}


