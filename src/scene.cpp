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
	
	shape.elapsedTime = elapsedTime;

	int const N_step = 4; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	int const N_stabilization = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	int const N_solver = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)


	for (int k_step = 0; simulation_running == true && k_step < N_step; ++k_step)
	{
		float dt_step = parameters.dt / N_step;
		simulation_compute_force(shape, parameters);

		cgp::buffer<vec3> prevX = {};
		for (const particle_element& particle : shape.particles) {
			prevX.push_back(particle.position);
		}
		// One step of numerical integration


		simulation_numerical_integration(shape, parameters, dt_step);
		for (int k_stabilization = 0; simulation_running == true && k_stabilization < N_stabilization; ++k_stabilization)
		{
			simulation_apply_constraints(shape, prevX, constraint, 1);
		}
		for (int k_solver = 0; simulation_running == true && k_solver < N_solver; ++k_solver)
		{
			preCalculations(shape);
			shapeMatching(shape, parameters, 1); //check parameters
		}
		adjustVelocity(shape, prevX, dt_step);

	}

	// Display particles
	if (gui.display_particles) {
		for (int k = 0; k < shape.particles.size(); ++k) {
			vec3 const& p = shape.particles[k].position;
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

	shape.initialize(0.3f, cgp::vec3(0.7, 1.3, 0.0), cgp::vec3(0, 75, 0));
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
		shape.initialize(0.3f, cgp::vec3(0.7, 1.3, 0.0), cgp::vec3(0,0,0));
		simulation_running = true;
	}

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