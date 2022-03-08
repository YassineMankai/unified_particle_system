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


// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt);

// Apply the constraints (fixed position, obstacles) on the cloth position and velocity
void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint);

// Helper function that tries to detect if the simulation diverged 
bool simulation_detect_divergence(cloth_structure const& cloth);


// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(shape_structure& shape, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(shape_structure& shape, simulation_parameters const& parameters, float dt);

void shapeMatching(shape_structure& shape, simulation_parameters const& parameters, float di);

void adjustVelocity(shape_structure& shape, cgp::buffer<cgp::vec3>& prevX, float dt);

// Apply the constraints (fixed position, obstacles) on the cloth position and velocity
void simulation_apply_constraints(shape_structure& shape, cgp::buffer<cgp::vec3>& prevX,  constraint_structure const& constraint, float di);

void calculateOptimalRotation(shape_structure& shape);

void calculateCurrentCom(shape_structure& shape);
void preCalculations(shape_structure& shape);