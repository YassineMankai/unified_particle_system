#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"
#include "../shape/custom_shape.hpp"
#include "../constraint/constraint.hpp"
#include "../uniform_grid/Rgrid.hpp"
#define EIGEN_NO_DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include "../third_party/eigen/Eigen/EigenValues"
#include "../third_party/eigen/Eigen/Eigen"
#include "custom_special_types.hpp"


struct simulation_parameters
{
    float dt = 0.014f;        // time step for the numerical integration
    float sphere_radius = 0.01f;        // damping parameter
    float alpha = 0.7;
    float beta = 0.5;
    int N_step = 12; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
    int N_stabilization = 2; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
    int N_solver = 2; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
    cgp::vec3 sphere3Pos = cgp::vec3(0.0f, -0.22f, 0.15f);
    float clothStiffness = 0.8f;
};

mat9 calculateInverseWithEigen(mat9 A);

// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(cgp::buffer<particle_element>& all_particles, float dt);


// Apply the constraints (fixed position, obstacles, shape matching, cloth constraints, fluid constraints, ...) 
void shapeMatching(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, float alpha, float beta);

void simulation_apply_shape_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes, constraint_structure const& constraint, simulation_parameters const& parameters);

void simulation_apply_contact_constraints(cgp::buffer<particle_element>& all_particles, const cgp::buffer<shape_structure>& all_shapes, std::vector<Rgrid> const& regularGrids, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float dt);
void simulation_apply_stabilization_contact_constraints(cgp::buffer<particle_element>& all_particles, const cgp::buffer<shape_structure>& all_shapes, std::vector<Rgrid> const& regularGrids, cgp::buffer<cgp::vec3>& prevX, constraint_structure const& constraint, float dt);

void calculateOptimalRotation(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);

void calculateCurrentCom(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);

void preCalculations(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);

// Calculate a position based velocity and update it with a collision correction term if needed
void adjustVelocity(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, float dt);