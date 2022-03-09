#include "scene.hpp"


using namespace cgp;

//void update_field_color(grid_2D<vec3>& field, buffer<particle_element> const& particles);

void scene_structure::display(double elapsedTime)
{
	// Basics common elements
	// ***************************************** //
	timer.update();
	environment.light = environment.camera.position();
	if (gui.display_frame)
		draw(global_frame, environment);



	// Elements of the scene: Obstacles (floor, sphere), and fixed position
	// ***************************************** //
	for (int wall_index = 0; wall_index < constraint.walls.size(); wall_index++) {
		const auto& wall = constraint.walls[wall_index];
		
		if (wall_index > 0) {
			obstacle_floor.transform.translation = 1.5* wall.normal;
			obstacle_floor.transform.rotation = rotation_transform::from_axis_angle(vec3(wall.normal.y, wall.normal.x,0), Pi / 2);
		}
		else {
			obstacle_floor.transform.translation = wall.point;
			obstacle_floor.transform.rotation = cgp::rotation_transform::from_axis_angle({ 0,1, 0},0);
		}
		draw(obstacle_floor, environment);
	}
	for (const auto& sphere : constraint.spheres)
	{
		obstacle_sphere.transform.translation = sphere.center;
		draw(obstacle_sphere, environment);
	}
	
	//shape.elapsedTime = elapsedTime;

	int const N_step = 4; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	int const N_stabilization = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	int const N_solver = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)


	for (int k_step = 0; simulation_running == true && k_step < N_step; ++k_step)
	{
		float dt_step = parameters.dt / N_step;
		simulation_compute_force(all_particles, parameters);

		cgp::buffer<vec3> prevX = {};
		for (const particle_element& particle : all_particles) {
			prevX.push_back(particle.position);
		}
		// One step of numerical integration


		simulation_numerical_integration(all_particles,dt_step);
		
		for (int k_stabilization = 0; simulation_running == true && k_stabilization < N_stabilization; ++k_stabilization)
		{
			simulation_apply_constraints(all_particles, prevX, constraint, 1);
		}
		
		for (int k_solver = 0; simulation_running == true && k_solver < N_solver; ++k_solver)
		{
			preCalculations(all_particles,all_shapes);
			shapeMatching(all_particles, all_shapes); //check parameters
		}
		
		adjustVelocity(all_particles, prevX, dt_step);
		

	}

	// Display particles
	if (gui.display_particles) {
		for (int k = 0; k < all_particles.size(); ++k) {
			vec3 const& p = all_particles[k].position;
			sphere_particle.transform.translation = p;
			draw(sphere_particle, environment);
		}
	}

	/*if (gui.display_color) {
		update_field_color(field, particles);
		opengl_update_texture_image(field_quad.texture, field);
		draw(field_quad, environment);
	}*/
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

	
	//cloth_texture = opengl_load_texture_image("assets/cloth.jpg");
	
	/*field.resize(30, 30);
	field_quad.initialize(mesh_primitive_quadrangle({ -1,-1,0 }, { 1,-1,0 }, { 1,1,0 }, { -1,1,0 }), "Field Quad");
	field_quad.shading.phong = { 1,0,0 };
	field_quad.texture = opengl_load_texture_image(field);*/

	//shape.initialize(0.3f, cgp::vec3(0.7, 1.3, 0.0), cgp::vec3(0, 75, 0));

	addCube(0.3f, cgp::vec3(0.7, 1.3, 0.0), cgp::vec3(0, 75, 0));

	addCube(0.3f, cgp::vec3(1.3, 1.3, 0.0), cgp::vec3(0, 75, 0));


	sphere_particle.initialize(mesh_primitive_sphere(), "Sphere particle");
	sphere_particle.transform.scaling = 0.01f;

}

void scene_structure::display_gui()
{
	bool reset = false;

	ImGui::Text("Display");
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Color", &gui.display_color);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Time step", &parameters.dt, 0.001f, 0.2f);
	ImGui::SliderFloat("Stiffness", &parameters.K, 1.0f, 80.0f, "%.3f", 2.0f);
	ImGui::SliderFloat("Wind magnitude", &parameters.wind.magnitude, 0, 60, "%.3f", 2.0f);
	ImGui::SliderFloat("Damping", &parameters.mu, 1.0f, 100.0f);
	ImGui::SliderFloat("Mass", &parameters.mass_total, 0.2f, 5.0f, "%.3f", 2.0f);

	ImGui::Spacing(); ImGui::Spacing();

	reset |= ImGui::SliderInt("Cloth samples", &gui.N_sample_edge, 4, 80);

	ImGui::Spacing(); ImGui::Spacing();
	reset |= ImGui::Button("Restart");
	if (reset) {
		all_particles.clear();
		all_shapes.clear();
		addCube(0.3f, cgp::vec3(0.7, 1.3, 0.0), cgp::vec3(0, 75, 0));
		addCube(0.3f, cgp::vec3(1.3, 1.3, 0.0), cgp::vec3(0, 75, 0));
		simulation_running = true;
	}

}

void scene_structure::addCube(float c, cgp::vec3 globalPosition, cgp::vec3 anglesEuler) {
	// Initial particle spacing (relative to h)
	float const h = p_parameters.h;
	// Fill a square with particles
	
	int indexShape = all_shapes.size();

	int indexStart = all_particles.size();



	shape_structure shape;

	shape.relativeLocationsOffset = indexStart;

	int sum = 0;
	cgp::vec3 centerOfMass = cgp::vec3(0.0, 0.0, 0.0);
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
				particle.phase = indexShape;
				particle.position = vec3(x_p + h / 8.0f, y_p + h / 8.0f, z_p + h / 8.0f) + vec3(0.0, 0.0, 1.0f) + globalPosition; // a zero value in z position will lead to a 2D simulation
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


	for (int i = indexStart; i < indexStart+sum; i++) {
		all_particles[i].position = (rotation * (all_particles[i].position - shape.com0)) + shape.com0;
	};


	/*
	vec3 res = vec3(0.0, 0.0, 0.0);
	for (int i = indexStart; i < indexStart + sum; i++) {
		particle_element& particle = all_particles[i];
		res = res + particle.position;
	};

	res = res / sum;

	if (norm(res - shape.com0) < 0.001) {
		for (int i = indexStart; i < indexStart + sum; i++) {
			particle_element& particle = all_particles[i];
			particle.position = particle.position - res + shape.com0;
		};
	}
	*/
	

	//initializing relative locations
	shape.relativeLocations = {};
	for (int i = indexStart; i <indexStart + sum; i++) {
		const particle_element& particleP = all_particles[i];
		shape.relativeLocations.push_back(particleP.position - shape.com0);
	}

	//pre calculation of Aqq
	cgp::mat3 p = cgp::mat3();
	cgp::mat3 q = cgp::mat3();
	for (int i = indexStart; i <indexStart + sum; i++) {
		const particle_element& particle = all_particles[i];
		p(0, 0) = shape.relativeLocations[i-indexStart].x;
		p(1, 0) = shape.relativeLocations[i-indexStart].y;
		p(2, 0) = shape.relativeLocations[i-indexStart].z;

		q(0, 0) = shape.relativeLocations[i-indexStart].x;
		q(0, 1) = shape.relativeLocations[i-indexStart].y;
		q(0, 2) = shape.relativeLocations[i-indexStart].z;

		shape.Aqq += particle.mass * p * q;
	}
	all_shapes.push_back(shape);
}

/*
void update_field_color(grid_2D<vec3>& field, buffer<particle_element> const& particles)
{
	field.fill({ 1,1,1 });
	float const d = 0.1f;
	int const Nf = int(field.dimension.x);
	for (int kx = 0; kx < Nf; ++kx) {
		for (int ky = 0; ky < Nf; ++ky) {

			float f = 0.0f;
			vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };
			for (size_t k = 0; k < particles.size(); ++k) {
				vec3 const& pi = particles[k].p;
				float const r = norm(pi - p0) / d;
				f += 0.25f * std::exp(-r * r);
			}
			field(kx, Nf - 1 - ky) = vec3(clamp(1 - f, 0, 1), clamp(1 - f, 0, 1), 1);
		}
	}
}
*/