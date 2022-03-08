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
	draw(obstacle_floor, environment);
	draw(obstacle_sphere, environment);
	for (auto const& c : constraint.fixed_sample)
	{
		sphere_fixed_position.transform.translation = c.second.position;
		draw(sphere_fixed_position, environment);
	}

	shape.elapsedTime = elapsedTime;
	// Simulation of the cloth
	// ***************************************** //
	int const N_step = 7; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	
	// iteration for constraints only
	
		// Update the forces on each particle
		//simulation_compute_force(cloth, parameters);

		// One step of numerical integration
		//simulation_numerical_integration(cloth, parameters, parameters.dt / N_step);

		// Apply the positional (and velocity) constraints
		//simulation_apply_constraints(cloth, constraint);




		simulation_compute_force(shape, parameters);
		
		cgp::buffer<vec3> prevX = {}; 
		for (const particle_element& particle : shape.particles) {
			prevX.push_back(particle.position);
		}
		// One step of numerical integration
		
		
		simulation_numerical_integration(shape, parameters, parameters.dt / N_step);
		//for (int k_step = 0; simulation_running == true && k_step < N_step; ++k_step)
		//{
			simulation_apply_constraints(shape, prevX, constraint);
		//}
		preCalculations(shape);
		shapeMatching(shape, parameters, parameters.dt / N_step); //check parameters

		adjustVelocity(shape, prevX, parameters.dt / N_step);
		

		/*// Check if the simulation has not diverged - otherwise stop it
		bool const simulation_diverged = simulation_detect_divergence(cloth);
		if (simulation_diverged) {
			std::cout << "\n *** Simulation has diverged ***" << std::endl;
			std::cout << " > The simulation is stoped" << std::endl;
			simulation_running = false;
		}*/
	


	// Cloth display
	// ***************************************** //

	// Prepare to display the updated cloth
	cloth.update_normal();        // compute the new normals
	cloth_drawable.update(cloth); // update the positions on the GPU

	// Display the cloth
	draw(cloth_drawable, environment);
	if (gui.display_wireframe)
		draw_wireframe(cloth_drawable, environment);


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

// Compute a new cloth in its initial position (can be called multiple times)
void scene_structure::initialize_cloth(int N_sample)
{
	cloth.initialize(N_sample);
	cloth_drawable.initialize(N_sample);
	cloth_drawable.drawable.texture = cloth_texture;


	constraint.fixed_sample.clear();
	constraint.add_fixed_position(0, 0, cloth);
	constraint.add_fixed_position(0, N_sample - 1, cloth);
}

void scene_structure::initialize()
{
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.look_at({ 3.0f,2.0f,2.0f }, { 0,0,0 }, { 0,0,1 });

	obstacle_floor.initialize(mesh_primitive_quadrangle({ -1.5f,-1.5f,0 }, { -1.5f,1.5f,0 }, { 1.5f,1.5f,0 }, { 1.5f,-1.5f,0 }));
	obstacle_floor.texture = opengl_load_texture_image("assets/wood.jpg");
	obstacle_floor.transform.translation = { 0,0,constraint.ground_z };

	obstacle_sphere.initialize(mesh_primitive_sphere());
	obstacle_sphere.transform.translation = constraint.sphere.center;
	obstacle_sphere.transform.scaling = constraint.sphere.radius;
	obstacle_sphere.shading.color = { 1,0,0 };

	sphere_fixed_position.initialize(mesh_primitive_sphere());
	sphere_fixed_position.transform.scaling = 0.02f;
	sphere_fixed_position.shading.color = { 0,0,1 };

	
	cloth_texture = opengl_load_texture_image("assets/cloth.jpg");
	initialize_cloth(gui.N_sample_edge);
	
	/*field.resize(30, 30);
	field_quad.initialize(mesh_primitive_quadrangle({ -1,-1,0 }, { 1,-1,0 }, { 1,1,0 }, { -1,1,0 }), "Field Quad");
	field_quad.shading.phong = { 1,0,0 };
	field_quad.texture = opengl_load_texture_image(field);*/

	shape.initialize(0.3f,cgp::vec3(0.7,1.3,0.0),75.f);
	sphere_particle.initialize(mesh_primitive_sphere(), "Sphere particle");
	sphere_particle.transform.scaling = 0.01f;

}

void scene_structure::display_gui()
{
	bool reset = false;

	ImGui::Text("Display");
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	ImGui::Checkbox("Texture Cloth", &cloth_drawable.drawable.shading.use_texture);
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Color", &gui.display_color);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Time step", &parameters.dt, 0.01f, 0.2f);
	ImGui::SliderFloat("Stiffness", &parameters.K, 1.0f, 80.0f, "%.3f", 2.0f);
	ImGui::SliderFloat("Wind magnitude", &parameters.wind.magnitude, 0, 60, "%.3f", 2.0f);
	ImGui::SliderFloat("Damping", &parameters.mu, 1.0f, 100.0f);
	ImGui::SliderFloat("Mass", &parameters.mass_total, 0.2f, 5.0f, "%.3f", 2.0f);

	ImGui::Spacing(); ImGui::Spacing();

	reset |= ImGui::SliderInt("Cloth samples", &gui.N_sample_edge, 4, 80);

	ImGui::Spacing(); ImGui::Spacing();
	reset |= ImGui::Button("Restart");
	if (reset) {
		initialize_cloth(gui.N_sample_edge);
		shape.initialize(0.3f, cgp::vec3(0.7, 1.3, 0.0), 75.f);
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