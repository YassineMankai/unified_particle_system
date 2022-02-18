#pragma once

#include "cgp/cgp.hpp"


// SPH Particle
struct particle_element
{
    cgp::vec3 position; // Position
    cgp::vec3 velocity; // Speed
    cgp::vec3 force; // Force

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

// Stores the buffers representing the cloth vertices
struct shape_structure
{    
    // Buffers are stored as 2D grid that can be accessed as grid(ku,kv)
    particle_parameters_structure p_parameters; //particle parameters
    cgp::buffer<particle_element> particles;      // Storage of the particles
    


    void initialize(float c = 0.3f);  // fill the shape with particles
};
