#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"
#include "../shape/shape.hpp"
#include "../constraint/constraint.hpp"
#define EIGEN_NO_DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include "../third_party/eigen/Eigen/EigenValues"
#include "../third_party/eigen/Eigen/Eigen"





struct simulation_parameters
{
    float dt = 0.01f;        // time step for the numerical integration
    float mass_total = 1.0f; // total mass of the cloth
    float K = 10.0f;         // stiffness parameter
    float mu = 40.0f;        // damping parameter

    //  Wind magnitude and direction
    struct {
        float magnitude = 0.0f;
        cgp::vec3 direction = { 0,-1,0 };
    } wind;
};


// Helper function that tries to detect if the simulation diverged 
bool simulation_detect_divergence(cloth_structure const& cloth);


// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(cgp::buffer<particle_element>& all_particles, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(cgp::buffer<particle_element>& all_particles, float dt);

void shapeMatching(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);

void adjustVelocity(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX, float dt);

// Apply the constraints (fixed position, obstacles) on the cloth position and velocity
void simulation_apply_constraints(cgp::buffer<particle_element>& all_particles, cgp::buffer<cgp::vec3>& prevX,  constraint_structure const& constraint, float di);

void calculateOptimalRotation(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);

void calculateCurrentCom(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);
void preCalculations(cgp::buffer<particle_element>& all_particles, cgp::buffer<shape_structure>& all_shapes);



