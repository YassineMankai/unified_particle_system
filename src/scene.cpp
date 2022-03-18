#include "scene.hpp"


using namespace cgp;

void scene_structure::setShapes(DemoScene demoType) {
	all_particles.clear();
	all_shapes.clear();
	constraint.fixed_sample.clear();

	switch (demoType) {
	case sc_LINEAR:
		parameters.dt = 0.018f;
		parameters.alpha = 0.8f;
		parameters.beta = 0.5f;
		parameters.N_step = 8;
		parameters.N_stabilization = 2;
		parameters.N_solver = 5;
		constraint.spheres = { {{-0.22f, 0.0f, 0.15f}, 0.15f},
								{ {0.22f, 0.0f, 0.15f}, 0.15f},
								{ {0.0f, -0.22f, 0.15f}, 0.15f},
								{ {0.0f, 0.22f, 0.15f}, 0.15f},
		};
		addCube(1.4f, 1.4f, 1.4f, 2, cgp::vec3(0.0f, 0.0f, 0.6f), cgp::vec3(0, 0, 0));
		break;
	case sc_QUADRATIC:
		parameters.dt = 0.020f;
		parameters.alpha = 0.8f;
		parameters.beta = 0.5f;
		parameters.N_step = 8;
		parameters.N_stabilization = 2;
		parameters.N_solver = 5;
		constraint.spheres = { {{-0.22f, 0.0f, 0.15f}, 0.15f},
								{ {0.22f, 0.0f, 0.15f}, 0.15f},
								{ {0.0f, -0.22f, 0.15f}, 0.15f},
								{ {0.0f, 0.22f, 0.15f}, 0.15f},
		};
		addCubeQuadratic(1.2f, 1.2f, 1.2f, 2, cgp::vec3(0.0f, 0.0f, 0.6f), cgp::vec3(0, 0, 0));
		break;
	case sc_CURTAIN:
		parameters.dt = 0.020f;
		parameters.alpha = 0.3f;
		parameters.beta = 0.15f;
		parameters.N_step = 8;
		parameters.N_stabilization = 2;
		parameters.N_solver = 4;
		parameters.clothStiffness = 0.9f;
		constraint.spheres = { {{-1.5f, -0.75f, 0.15f}, 0.15f},
								{ {-1.5f, -0.25f, 0.15f}, 0.15f},
								{ {0.2f, 0.0f, 0.4f}, 0.15f},
								{ {-1.5f, 0.75f, 0.15f}, 0.15f},
		};
		addCube(1.0f, 1.0f, 1.0f, 2, cgp::vec3(0.05f, 0.0f, 0.7f), cgp::vec3(0, Pi / 4, 0));
		addCloth(3.0f, 4.0f, 3, cgp::vec3(-0.05f, 0.0f, 0.5f), cgp::vec3(0, 0, 0), 1);
		break;
	case sc_TRAMPOLINE:
		parameters.dt = 0.02f;
		parameters.alpha = 0.8f;
		parameters.beta = 0.5f;
		parameters.N_step = 10;
		parameters.N_stabilization = 4;
		parameters.N_solver = 4;
		parameters.clothStiffness = 0.9f;
		constraint.spheres = { {{-1.5f, -0.75f, 0.15f}, 0.15f},
								{ {-1.5f, -0.25f, 0.15f}, 0.15f},
								{ {0.2f, 0.0f, 0.3f}, 0.15f},
								{ {-1.5f, 0.75f, 0.15f}, 0.15f},
		};
		addCube(0.8f, 0.8f, 0.8f, 2, cgp::vec3(-0.25f, 0.0f, 1.2f), cgp::vec3(Pi / 4, Pi / 4, 0));
		addCloth(3.0f, 3.0f, 2, cgp::vec3(-0.2f, 0.0f, 0.4f), cgp::vec3(0, Pi / 2 + Pi / 12, 0), 2);
		break;
	}
	parameters.sphere3Pos = constraint.spheres[2].center;
}

void scene_structure::initialize()
{
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.look_at({ 3.0f,2.0f,2.0f }, { 0,0,0 }, { 0,0,1 });

	obstacle_floor.initialize(mesh_primitive_quadrangle({ -1.5f,-1.5f,0 }, { -1.5f,1.5f,0 }, { 1.5f,1.5f,0 }, { 1.5f,-1.5f,0 }));
	obstacle_floor.texture = opengl_load_texture_image("assets/wood.jpg");
	obstacle_floor.transform.translation = { 0,0,0 }; //TODO display floor

	obstacle_sphere.initialize(mesh_primitive_sphere());
	obstacle_sphere.transform.translation = { 0.1f, 0.5f, 0.0f };
	obstacle_sphere.transform.scaling = 0.15f;
	obstacle_sphere.shading.color = { 1,0,0 };

	setShapes();


	sphere_particle.initialize(mesh_primitive_sphere(), "Sphere particle");
	sphere_particle.transform.scaling = parameters.sphere_radius;
}

void scene_structure::simulate() {
	for (int k_step = 0; simulation_running == true && k_step < parameters.N_step; ++k_step)
	{
		// ****************************** //
		// computing forces
		// ****************************** //
		float dt_step = parameters.dt / parameters.N_step;
		simulation_compute_force(all_particles, all_shapes, parameters);

		// ****************************** //
		// storing previous positions
		// ****************************** //
		cgp::buffer<vec3> prevX = {};
		for (const particle_element& particle : all_particles) {
			prevX.push_back(particle.position);
		}
		// One step of numerical integration

		// ****************************** //
		// predict positions
		// ****************************** //
		simulation_numerical_integration(all_particles, dt_step);

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

		// ****************************** //
		// Solve contact constraints
		// ****************************** //

		for (int k_stabilization = 0; simulation_running == true && k_stabilization < parameters.N_stabilization; ++k_stabilization)
		{
			simulation_apply_contact_constraints(all_particles, all_shapes, regularGrids, prevX, constraint, dt_step / parameters.N_stabilization);
		}

		// ****************************** //
		// Solve shape constraints
		// ****************************** //

		for (int k_solver = 0; simulation_running == true && k_solver < parameters.N_solver; ++k_solver)
		{
			preCalculations(all_particles, all_shapes);
			shapeMatching(all_particles, all_shapes, parameters.alpha, parameters.beta); //check parameters
			simulation_apply_shape_constraints(all_particles, all_shapes, constraint, parameters);
			simulation_apply_contact_constraints(all_particles, all_shapes, regularGrids, prevX, constraint, dt_step / parameters.N_solver);
		}

		//put here the constraints of the cloth for example

		// ****************************** //
		// Coorecting velocity
		// ****************************** //

		adjustVelocity(all_particles, prevX, dt_step);

	}
}

void scene_structure::display(double elapsedTime)
{
	// Basics common elements
	// ***************************************** //
	timer.update();
	environment.light = environment.camera.position();
	if (gui.display_frame)
		draw(global_frame, environment);


	// Elements of the scene: Obstacles (floor, sphere), and fixed position
	// * //
	for (int wall_index = 0; wall_index < constraint.walls.size(); wall_index++) {
		const auto& wall = constraint.walls[wall_index];

		if (wall_index > 0) {
			obstacle_floor.transform.translation = 1.5 * wall.normal - cgp::vec3(0, 0, 1);
			obstacle_floor.transform.rotation = rotation_transform::from_axis_angle(vec3(wall.normal.y, wall.normal.x, 0), Pi / 2);
		}

		else {
			obstacle_floor.transform.translation = wall.point;
			obstacle_floor.transform.rotation = cgp::rotation_transform::from_axis_angle({ 0,1, 0 }, 0);
		}

		draw(obstacle_floor, environment);
	}

	constraint.spheres[2].center = parameters.sphere3Pos;
	for (const auto& sphere : constraint.spheres)
	{
		obstacle_sphere.transform.translation = sphere.center;
		obstacle_sphere.transform.scaling = sphere.radius;
		draw(obstacle_sphere, environment);
	}

	simulate();

	// Display particles
	for (int k = 0; k < all_particles.size(); ++k) {
		vec3 const& p = all_particles[k].position;
		sphere_particle.transform.translation = p;
		draw(sphere_particle, environment);
	}

}

void scene_structure::display_gui()
{

	ImGui::Text("Display");
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Time step", &parameters.dt, 0.001f, 0.2f);
	ImGui::SliderInt("N_step", &parameters.N_step, 1, 12, "%100.0f");
	ImGui::SliderInt("N_stabilization", &parameters.N_stabilization, 1, 8, "%100.0f");
	ImGui::SliderInt("N_solver", &parameters.N_solver, 1, 8, "%100.0f");

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::SliderFloat("alpha", &parameters.alpha, 0.05, 1.5, "%.3g", 1);
	ImGui::SliderFloat("beta", &parameters.beta, 0.05, 1, "%.3g", 1);
	ImGui::Spacing();
	ImGui::SliderFloat("cloth stiffness", &parameters.clothStiffness, 0.01, 1.5f, "%.3f", 1.0f);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::SliderFloat("SphereX", &parameters.sphere3Pos.x, -1.8, 1.8, "%.05f", 1.0f);
	ImGui::SliderFloat("SphereY", &parameters.sphere3Pos.y, -1.8, 1.8, "%.05f", 1.0f);
	ImGui::SliderFloat("SphereZ", &parameters.sphere3Pos.z, -0.2, 1.8, "%.05f", 1.0f);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Spacing(); ImGui::Spacing();



	bool resetWithLinear = false;
	resetWithLinear |= ImGui::Button("Reset with Linear");
	if (resetWithLinear) {
		setShapes(sc_LINEAR);
		simulation_running = false;
	}

	bool resetWithQuadratic = false;
	resetWithQuadratic |= ImGui::Button("Reset with Quadratic");
	if (resetWithQuadratic) {
		setShapes(sc_QUADRATIC);
		simulation_running = false;
	}

	bool resetWithCurtain = false;
	resetWithCurtain |= ImGui::Button("Reset with Curtain");
	if (resetWithCurtain) {
		setShapes(sc_CURTAIN);
		simulation_running = false;
	}

	bool resetWithTrampoline = false;
	resetWithTrampoline |= ImGui::Button("Reset with Trampoline");
	if (resetWithTrampoline) {
		setShapes(sc_TRAMPOLINE);
		simulation_running = false;
	}

	if (ImGui::Button("Start/Stop"))
		simulation_running = !simulation_running;
}

void scene_structure::addCube(float c_x, float c_y, float c_z, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler) {
	// Initial particle spacing (relative to h)
	float const h = scene_parameters.shape_size;
	// Fill a square with particles

	int indexShape = all_shapes.size();
	int indexStart = all_particles.size();

	shape_structure shape;
	shape.type = RIGID;

	shape.relativeLocationsOffset = indexStart;

	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);

	for (float x = -h * c_x; x <= h * c_x; x = x + (spacing * 0.01f))
	{
		for (float y = -h * c_y; y <= h * c_y; y = y + (spacing * 0.01f))
		{
			for (float z = -h * c_z; z <= h * c_z; z = z + (spacing * 0.01f))
			{
				float x_p = x;
				float y_p = y;
				float z_p = z;
				particle_element particle;
				particle.phase = indexShape;
				particle.position = vec3(x_p, y_p, z_p); // a zero value in z position will lead to a 2D simulation
				centerOfMass += particle.position;
				sum++;
				all_particles.push_back(particle);
			}
		}
	}

	//initializing center of mass
	centerOfMass /= sum;
	shape.com0 = centerOfMass;
	shape.nbParticles = sum;

	//initializing relative locations
	shape.relativeLocations = {};
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particleP = all_particles[i];
		shape.relativeLocations.push_back(particleP.position - shape.com0);
	}

	//pre calculation of Aqq
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particle = all_particles[i];
		p(0, 0) = shape.relativeLocations[i - indexStart].x;
		p(1, 0) = shape.relativeLocations[i - indexStart].y;
		p(2, 0) = shape.relativeLocations[i - indexStart].z;

		q(0, 0) = shape.relativeLocations[i - indexStart].x;
		q(0, 1) = shape.relativeLocations[i - indexStart].y;
		q(0, 2) = shape.relativeLocations[i - indexStart].z;

		shape.Aqq += particle.mass * p * q;
	}
	shape.Aqq = inverse(shape.Aqq);

	rotation_transform rotationy = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), anglesEuler.y);
	rotation_transform rotationx = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), anglesEuler.x);
	rotation_transform rotationz = rotation_transform::from_axis_angle(vec3(0.0, 0.0, 1.0), anglesEuler.z);
	rotation_transform rotation = rotationx * rotationy * rotationz;
	affine_rt rt = rotation_around_center(rotation, shape.com0);
	for (int i = indexStart; i < indexStart + sum; i++) {
		all_particles[i].position = rt * all_particles[i].position + globalPosition;
		//if (i - indexStart < 3)
			//constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(i, all_particles[i].position + vec3(rand_interval(0, 0.03), rand_interval(0, 0.03), rand_interval(0, 0.03))));
	};

	all_shapes.push_back(shape);
}
void scene_structure::addCubeQuadratic(float c_x, float c_y, float c_z, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler)
{
	// Initial particle spacing (relative to h)
	float const h = scene_parameters.shape_size;
	// Fill a square with particles

	int indexShape = all_shapes.size();

	int indexStart = all_particles.size();
	shape_structure shape;
	shape.type = QUADRATIC;
	shape.relativeLocationsOffset = indexStart;
	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);
	for (float x = -h * c_x; x <= h * c_x; x = x + (spacing * 0.01f))
	{
		for (float y = -h * c_y; y <= h * c_y; y = y + (spacing * 0.01f))
		{
			for (float z = -h * c_z; z <= h * c_z; z = z + (spacing * 0.01f))
			{
				float x_p = x;
				float y_p = y;
				float z_p = z;
				particle_element particle;
				particle.phase = indexShape;
				particle.position = vec3(x_p, y_p, z_p); // a zero value in z position will lead to a 2D simulation
				centerOfMass += particle.position;
				sum++;
				all_particles.push_back(particle);
			}
		}
	}


	//initializing center of mass
	centerOfMass /= sum;
	shape.com0 = centerOfMass;
	shape.nbParticles = sum;

	//initializing relative locations
	shape.relativeLocations = {};
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particleP = all_particles[i];
		shape.relativeLocations.push_back(particleP.position - shape.com0);
	}

	//pre calculation of Aqq
	cgp::mat9 p = cgp::mat9();
	cgp::mat9 q = cgp::mat9();

	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particle = all_particles[i];
		p(0, 0) = shape.relativeLocations[i - indexStart].x;
		p(1, 0) = shape.relativeLocations[i - indexStart].y;
		p(2, 0) = shape.relativeLocations[i - indexStart].z;

		p(3, 0) = pow(shape.relativeLocations[i - indexStart].x, 2);
		p(4, 0) = pow(shape.relativeLocations[i - indexStart].y, 2);
		p(5, 0) = pow(shape.relativeLocations[i - indexStart].z, 2);

		p(6, 0) = shape.relativeLocations[i - indexStart].x * shape.relativeLocations[i - indexStart].y;
		p(7, 0) = shape.relativeLocations[i - indexStart].y * shape.relativeLocations[i - indexStart].z;
		p(8, 0) = shape.relativeLocations[i - indexStart].z * shape.relativeLocations[i - indexStart].x;

		cgp::vec9 v = cgp::vec9();

		v(0) = p(0, 0);
		v(1) = p(1, 0);
		v(2) = p(2, 0);

		v(3) = p(3, 0);
		v(4) = p(4, 0);
		v(5) = p(5, 0);

		v(6) = p(6, 0);
		v(7) = p(7, 0);
		v(8) = p(8, 0);

		shape.qQuad.push_back(v);

		q(0, 0) = shape.relativeLocations[i - indexStart].x;
		q(0, 1) = shape.relativeLocations[i - indexStart].y;
		q(0, 2) = shape.relativeLocations[i - indexStart].z;

		q(0, 3) = pow(shape.relativeLocations[i - indexStart].x, 2);
		q(0, 4) = pow(shape.relativeLocations[i - indexStart].y, 2);
		q(0, 5) = pow(shape.relativeLocations[i - indexStart].z, 2);

		q(0, 6) = shape.relativeLocations[i - indexStart].x * shape.relativeLocations[i - indexStart].y;
		q(0, 7) = shape.relativeLocations[i - indexStart].y * shape.relativeLocations[i - indexStart].z;
		q(0, 8) = shape.relativeLocations[i - indexStart].z * shape.relativeLocations[i - indexStart].x;

		shape.AqqQuad += particle.mass * p * q;
	}

	shape.AqqQuad = calculateInverseWithEigen(shape.AqqQuad);

	rotation_transform rotationy = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), anglesEuler.y);
	rotation_transform rotationx = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), anglesEuler.x);
	rotation_transform rotationz = rotation_transform::from_axis_angle(vec3(0.0, 0.0, 1.0), anglesEuler.z);
	rotation_transform rotation = rotationx * rotationy * rotationz;
	affine_rt rt = rotation_around_center(rotation, shape.com0);


	for (int i = indexStart; i < indexStart + sum; i++) {
		all_particles[i].position = rt * all_particles[i].position + globalPosition;
		//if (i - indexStart < 3)
			//constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(i, all_particles[i].position + vec3(rand_interval(0, 0.03), rand_interval(0, 0.03), rand_interval(0, 0.03))));
	};

	all_shapes.push_back(shape);
}
void scene_structure::addSphere(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler) {
	// Initial particle spacing (relative to h)
	float const h = scene_parameters.shape_size;
	// Fill a square with particles

	int indexShape = all_shapes.size();

	int indexStart = all_particles.size();



	shape_structure shape;
	shape.type = RIGID;

	shape.relativeLocationsOffset = indexStart;

	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);
	for (float x = -2 * h * c; x <= 2 * h * c; x = x + 0.45f * h)
	{
		for (float y = -2 * h * c; y <= 2 * h * c; y = y + 0.45f * h)
		{
			for (float z = -2 * h * c; z <= 2 * h * c; z = z + 0.45f * h)
			{
				float x_p = x;
				float y_p = y;
				float z_p = z;

				if (std::sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) < 1.25 * h) {
					particle_element particle;
					particle.phase = indexShape;
					particle.position = vec3(x_p, y_p, z_p) + globalPosition; // a zero value in z position will lead to a 2D simulation
					centerOfMass += particle.position;
					sum++;
					all_particles.push_back(particle);
				}
			}
		}
	}
	//initializing center of mass
	centerOfMass /= sum;
	shape.com0 = centerOfMass;
	shape.nbParticles = sum;

	rotation_transform rotationy = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), anglesEuler.y);
	rotation_transform rotationx = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), anglesEuler.x);
	rotation_transform rotationz = rotation_transform::from_axis_angle(vec3(0.0, 0.0, 1.0), anglesEuler.z);

	rotation_transform rotation = rotationx * rotationy * rotationz;

	affine_rt rt = rotation_around_center(rotation, shape.com0);


	for (int i = indexStart; i < indexStart + sum; i++) {
		all_particles[i].position = (rotation * (all_particles[i].position - shape.com0)) + shape.com0;
	};


	//initializing relative locations
	shape.relativeLocations = {};
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particleP = all_particles[i];
		shape.relativeLocations.push_back(particleP.position - shape.com0);
	}

	//pre calculation of Aqq
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particle = all_particles[i];
		p(0, 0) = shape.relativeLocations[i - indexStart].x;
		p(1, 0) = shape.relativeLocations[i - indexStart].y;
		p(2, 0) = shape.relativeLocations[i - indexStart].z;

		q(0, 0) = shape.relativeLocations[i - indexStart].x;
		q(0, 1) = shape.relativeLocations[i - indexStart].y;
		q(0, 2) = shape.relativeLocations[i - indexStart].z;

		shape.Aqq += particle.mass * p * q;
	}
	shape.Aqq = inverse(shape.Aqq);
	all_shapes.push_back(shape);
}
void scene_structure::addPyramid(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler) {
	// Initial particle spacing (relative to h)
	float const h = scene_parameters.shape_size * c;
	// Fill a square with particles
	int indexShape = all_shapes.size();
	int indexStart = all_particles.size();
	shape_structure shape;
	shape.type = RIGID;
	shape.relativeLocationsOffset = indexStart;

	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);
	for (float z = 0; z <= h; z = z + (2 * 0.01f))
	{
		for (float x = -h + z; x <= h - z; x = x + (2 * 0.01f))
		{
			for (float y = -h + z; y <= h - z; y = y + (2 * 0.01f))
			{
				float x_p = x;
				float y_p = y;
				float z_p = z;
				particle_element particle;
				particle.phase = indexShape;
				particle.position = vec3(x_p, y_p, z_p) + globalPosition; // a zero value in z position will lead to a 2D simulation
				centerOfMass += particle.position;
				sum++;
				all_particles.push_back(particle);
			}
		}
	}
	//initializing center of mass
	centerOfMass /= sum;
	shape.com0 = centerOfMass;
	shape.nbParticles = sum;

	rotation_transform rotationy = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), anglesEuler.y);
	rotation_transform rotationx = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), anglesEuler.x);
	rotation_transform rotationz = rotation_transform::from_axis_angle(vec3(0.0, 0.0, 1.0), anglesEuler.z);

	rotation_transform rotation = rotationx * rotationy * rotationz;

	affine_rt rt = rotation_around_center(rotation, shape.com0);


	for (int i = indexStart; i < indexStart + sum; i++) {
		all_particles[i].position = (rotation * (all_particles[i].position - shape.com0)) + shape.com0;
	};


	//initializing relative locations
	shape.relativeLocations = {};
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particleP = all_particles[i];
		shape.relativeLocations.push_back(particleP.position - shape.com0);
	}

	//pre calculation of Aqq
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	for (int i = indexStart; i < indexStart + sum; i++) {
		const particle_element& particle = all_particles[i];
		p(0, 0) = shape.relativeLocations[i - indexStart].x;
		p(1, 0) = shape.relativeLocations[i - indexStart].y;
		p(2, 0) = shape.relativeLocations[i - indexStart].z;

		q(0, 0) = shape.relativeLocations[i - indexStart].x;
		q(0, 1) = shape.relativeLocations[i - indexStart].y;
		q(0, 2) = shape.relativeLocations[i - indexStart].z;

		shape.Aqq += particle.mass * p * q;
	}
	shape.Aqq = inverse(shape.Aqq);
	all_shapes.push_back(shape);
}
void scene_structure::addCloth(float c_w, float c_h, int spacing, cgp::vec3 globalPosition, cgp::vec3 anglesEuler, int type) {
	// Initial particle spacing (relative to h)
	float const h = scene_parameters.shape_size;
	// Fill a square with particles
	int indexShape = all_shapes.size();
	int indexStart = all_particles.size();
	shape_structure shape;
	shape.relativeLocationsOffset = indexStart;
	shape.type = CLOTH;
	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);

	float step = (spacing * 0.01f);

	for (float z = -h * c_h; z <= h * c_h; z = z + step)
	{
		for (float y = -h * c_w; y <= h * c_w; y = y + step)
		{
			float y_p = y;
			float z_p = z;
			particle_element particle;
			particle.phase = indexShape;
			particle.position = vec3(0.0, y_p, z_p) + globalPosition; // a zero value in z position will lead to a 2D simulation
			centerOfMass += particle.position;
			sum++;
			all_particles.push_back(particle);
		}
	}

	for (float z = -h * c_h; z <= h * c_h; z = z + step)
	{
		shape.height++;
	}

	for (float y = -h * c_w; y <= h * c_w; y = y + step) {
		shape.width++;
	}

	//initializing center of mass
	centerOfMass /= sum;
	shape.com0 = centerOfMass;
	shape.nbParticles = sum;

	rotation_transform rotationy = rotation_transform::from_axis_angle(vec3(0.0, 1.0, 0.0), anglesEuler.y);
	rotation_transform rotationx = rotation_transform::from_axis_angle(vec3(1.0, 0.0, 0.0), anglesEuler.x);
	rotation_transform rotationz = rotation_transform::from_axis_angle(vec3(0.0, 0.0, 1.0), anglesEuler.z);
	rotation_transform rotation = rotationx * rotationy * rotationz;
	affine_rt rt = rotation_around_center(rotation, shape.com0);


	for (int i = indexStart; i < indexStart + sum; i++) {
		all_particles[i].position = (rotation * (all_particles[i].position - shape.com0)) + shape.com0;
	};

	////inserting fixed points
	if (type == 0) {
		int index1 = shape.getGlobalIndex(shape.height - 1, 0);
		int index2 = shape.getGlobalIndex(shape.height - 1, shape.width - 1);
		const particle_element& p1 = all_particles[index1];
		const particle_element& p2 = all_particles[index2];
		constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index1, p1.position));
		constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index2, p2.position));
	}
	else if (type == 2) {
		for (int i = 0; i < shape.width; i++) {

			int index1 = shape.getGlobalIndex(shape.height - 1, i);
			int index2 = shape.getGlobalIndex(0, i);

			const particle_element& p1 = all_particles[index1];
			const particle_element& p2 = all_particles[index2];

			constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index1, p1.position));
			constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index2, p2.position));

		}
		for (int i = 0; i < shape.height; i++) {

			int index1 = shape.getGlobalIndex(i, shape.width - 1);
			int index2 = shape.getGlobalIndex(i, 0);

			const particle_element& p1 = all_particles[index1];
			const particle_element& p2 = all_particles[index2];

			constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index1, p1.position));
			constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index2, p2.position));

		}
	}
	else {
		for (int i = 0; i < shape.width; i++) {

			int index1 = shape.getGlobalIndex(shape.height - 1, i);
			const particle_element& p1 = all_particles[index1];
			constraint.fixed_sample.insert(std::pair<int, cgp::vec3>(index1, p1.position));
		}
	}

	//initializing lengths of springs 
	int j = shape.width / 2;
	int i = shape.height / 2;
	int globalI = shape.getGlobalIndex(i, j);

	shape.structLength0 = norm(all_particles[globalI].position - all_particles[shape.getNeighbour(globalI, 0, 1)].position);
	shape.shearLength0 = norm(all_particles[globalI].position - all_particles[shape.getNeighbour(globalI, 1, 1)].position);
	shape.bendLength0 = norm(all_particles[globalI].position - all_particles[shape.getNeighbour(globalI, 2, 0)].position);

	all_shapes.push_back(shape);
}


