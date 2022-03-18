#include "simulation.hpp"

using namespace cgp;







// ################# Constraints #################

// Helper functions
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
			
			for (int j = 0; j < 3; j++) {
				p(j, 0) = particle.position[j] - shape.com[j];
				q(0,j) = shape.relativeLocations[i - shape.relativeLocationsOffset][j];
			}


			shape.Apq += particle.mass * p * q;
		}
		else if (shape.type == QUADRATIC) {
			
			for (int j = 0; j < 3; j++) {
				p(j, 0) = particle.position[j] - shape.com[j];
				q(0, j) = shape.relativeLocations[i - shape.relativeLocationsOffset][j];
			}

			shape.Apq += particle.mass * p * q;

			cgp::vec9& q = shape.qQuad[i - shape.relativeLocationsOffset];

			for (int j = 0; j < 3; j++) {
				pQuad(j, 0) = particle.position[j] - shape.com[j];
			}


			for (int i = 0; i < 9; i++) {
				qQuad(0, i) = q(i);
			}


			shape.ApqQuad += particle.mass * pQuad * qQuad;
		}
	}

	for (int i = 0; i < all_shapes.size(); i++) {
		shape_structure& shape = all_shapes[i];
		if (shape.type == RIGID) {
			cgp::mat3 sym = transpose(shape.Apq) * shape.Apq;
			Eigen::MatrixXd A(3, 3);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					A(i, j) = sym(i, j);
				}
			}
			Eigen::EigenSolver <Eigen::MatrixXd > es(A);
			Eigen::MatrixXcd Dc = es.eigenvalues().asDiagonal();
			Eigen::MatrixXcd Vc = es.eigenvectors();
			Eigen::MatrixXd D = Dc.real();
			Eigen::MatrixXd V = Vc.real();

			for (int j = 0; j < 3; j++) {
				D(j,j) = 1.f / sqrt(D(j, j));
			}

			Eigen::MatrixXd Sinv = (V * D * V.inverse());
			cgp::mat3 S_ = cgp::mat3();

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					S_(i, j) = Sinv(i, j);
				}
			}
			cgp::mat3 result = shape.Apq * S_;
			shape.optimalRotation = result;
			shape.A = shape.Apq * shape.Aqq;
			shape.A = shape.A / pow(cgp::det(shape.A), 1.f / 3.f);
		}
		else if (shape.type == QUADRATIC) {
			cgp::mat3 sym = transpose(shape.Apq) * shape.Apq;
			Eigen::MatrixXd A(3, 3);

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					A(i, j) = sym(i, j);
				}
			}
			Eigen::EigenSolver <Eigen::MatrixXd > es(A);
			Eigen::MatrixXcd Dc = es.eigenvalues().asDiagonal();
			Eigen::MatrixXcd Vc = es.eigenvectors();
			Eigen::MatrixXd D = Dc.real();
			Eigen::MatrixXd V = Vc.real();
			for (int j = 0; j < 3; j++) {
				D(j, j) = 1.f / sqrt(D(j, j));
			}
			
			Eigen::MatrixXd Sinv = (V * D * V.inverse());
			cgp::mat3 S_ = cgp::mat3();
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					S_(i, j) = Sinv(i, j);
				}
			}
			cgp::mat3 result = shape.Apq * S_;

			cgp::mat39& rotation = shape.optimalRotationQuad;

			rotation = cgp::mat39();
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					rotation(i, j) = result(i, j);
				}
			}
			cgp::mat39 aQuad = shape.ApqQuad * shape.AqqQuad;

			cgp::mat3 A_ = cgp::mat3();

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					A_(i, j) = aQuad(i, j);
				}
			}
			A_ = A_ / pow(cgp::det(A_), 1.f / 3.f);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					aQuad(i, j) = A_(i, j);
				}
			}
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

			cgp::buffer<float> Ls = { shape.structLength0,shape.structLength0,shape.structLength0,shape.structLength0,
									shape.shearLength0,shape.shearLength0 ,shape.shearLength0,shape.shearLength0,
									shape.bendLength0 ,shape.bendLength0 ,shape.bendLength0 ,shape.bendLength0 };
			cgp::buffer<int> is = {
				1,-1,0,0
				,1,-1,1,-1
				,0,0,2,-2 };

			cgp::buffer<int> js = {
				0,0,1,-1,
				1,-1,-1,1,
				2,-2,0,0 };


			for (int h = 0; h < Ls.size(); h++) {
				int index = shape.getNeighbour(pIndex, is[h], js[h]);

				if (index != -1) {
					const particle_element& p = all_particles[index];
					float  distance = norm(p.position - particle.position);
					dx += k * (distance - Ls[h]) * normalize(p.position - particle.position);
					numberOfConstraints++;
				}
			}
			x += dx / numberOfConstraints;

		}
	}
}

// Constraints related to the environment (walls, obstacles, fixed points, particle-particle)
void simulation_apply_contact_constraints(cgp::buffer<particle_element>& all_particles, const cgp::buffer<shape_structure>& all_shapes, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float dt)
{
	// update neighborhood data
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

	// particle-env contacts
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
				all_particles[pIndex].dx_friction_and_restitution += (-0.7 * vn + 0.4 * vpar - v) * dt; 
				all_particles[pIndex].nbConstraintFrictionRestitution += 1;
			}
		}
		for (const auto& sphere : constraint.spheres) {

			if (norm(sphere.center - p) < sphere.radius + 0.01) {
				vec3 dir = (p - sphere.center) / norm(p - sphere.center);
				dx += sphere.center + (sphere.radius + 0.01) * dir - p;
				all_particles[pIndex].nbConstraint += 1;
				vec3 vn = dot(v, dir) * dir;
				vec3 vpar = v - vn;
				all_particles[pIndex].dx_friction_and_restitution += (-0.7 * vn + 0.4 * vpar - v) * dt; 
				all_particles[pIndex].nbConstraintFrictionRestitution += 1;
			}
		}
	}

	// particle-particle contact
	for (int p1 = 0; p1 < all_particles.size(); p1++)
	{
		particle_element& particle1 = all_particles[p1];
		
		// retrieve neighbor particles
		std::vector<int> neighborhood;
		for (int i = 0; i < all_shapes.size(); i++) {
			std::vector<int> res = regularGrids[i].getNeighborhood(particle1.position);
			neighborhood.insert(neighborhood.end(), res.begin(), res.end());
		}

		// project collision constraint
		float const epsilon = 1e-5f;
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
				particle2.dx += -dx;

				particle1.nbConstraint += 1;
				particle2.nbConstraint += 1;
				
				particle1.nbConstraintFrictionRestitution += 1;
				particle2.nbConstraintFrictionRestitution += 1;

				if (norm(dx) > 0.01 * dt) {
					particle1.dx_friction_and_restitution += dx;
					particle2.dx_friction_and_restitution += -dx;
				}
				else
				{
					particle1.dx_friction_and_restitution += -particle1.velocity * dt / 3.f;
					particle2.dx_friction_and_restitution += -particle2.velocity * dt / 3.f;
				}
			}
		}
	}

	// Apply dx
	for (int i = 0; i < all_particles.size(); i++) {
		particle_element& particle = all_particles[i];
		if (particle.nbConstraint > 0) {
			particle.position += particle.dx / particle.nbConstraint;
			prevX[i] += particle.dx / particle.nbConstraint;
			if (particle.nbConstraintFrictionRestitution > 0) {
				prevX[i] -= particle.dx_friction_and_restitution / particle.nbConstraintFrictionRestitution;
				particle.nbConstraintFrictionRestitution = 0;
				particle.dx_friction_and_restitution = vec3(0, 0, 0);
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
	for (int pIndex = 0; pIndex < N; ++pIndex) {
		all_particles[pIndex].force = all_particles[pIndex].mass * g;
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


