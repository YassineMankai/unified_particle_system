#pragma once

#include "cgp/cgp.hpp"


// SPH Particle
struct particle_element
{
    cgp::vec3 position; // Position
    cgp::vec3 velocity; // Speed
    cgp::vec3 force; // Force
    bool flag = false;
    cgp::vec3 dv; // Force
    float mass = 0.0001;
    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : position{ 0,0,0 }, velocity{ 0,0,0 }, force{ 0,0,0 }, rho(0), pressure(0) {}
};

// SPH simulation parameters
struct particle_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.10f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

    // Total mass of a particle (consider rho0 h^2)
    float m = rho0 * h * h;

    // viscosity parameter
    float nu = 0.2f;

    // Stiffness converting density to pressure
    float stiffness = 8.0f;

};

// Stores the buffers representing the shape vertices
struct shape_structure
{    
    double elapsedTime = 0.0;
    // Buffers are stored as 2D grid that can be accessed as grid(ku,kv)
    particle_parameters_structure p_parameters; //particle parameters
    cgp::buffer<particle_element> particles;      // Storage of the particles
    cgp::vec3 com0; //intial center of mass
    cgp::buffer<cgp::vec3> relativeLocations;
    cgp::mat3 optimalRotation; //current optimal rotation
    cgp::mat3 A;//optimal linear transformation
    cgp::mat3 Aqq;
    cgp::vec3 com; //current center of mass
    void initialize(float c = 0.3f, cgp::vec3 globalPosition = cgp::vec3(0.0, 0.0, 0.0), float angle=0.f);  // fill the shape with particles
};
