#pragma once

#include "cgp/cgp.hpp"



// SPH Particle
struct particle_element
{
    cgp::vec3 position; // Position
    cgp::vec3 velocity; // Speed
    cgp::vec3 force; // Force
    bool flagConstraint = false;
    cgp::vec3 dv; // Velocity adjustment due to contact collision
    float mass = 0.0001;
    int phase;
    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : position{ 0,0,0 }, velocity{ 0,0,0 }, force{ 0,0,0 }, rho(0), pressure(0) {}
};

// SPH simulation parameters
struct particle_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.06f;

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
    int relativeLocationsOffset;
    int nbParticles;
    bool isQuadratic = false;
    cgp::buffer<cgp::vec3> relativeLocations;
    cgp::vec3 com; //current center of mass
    /// <summary>
    /// linear variables
    /// </summary>
    cgp::mat3 optimalRotation; //current optimal rotation
    cgp::mat3 A;//optimal linear transformation
    cgp::mat3 Aqq=cgp::mat3();
    cgp::mat3 Apq=cgp::mat3();
    

    //quadratic variables

    cgp::mat39 AQuad;
    cgp::mat39 optimalRotationQuad;
    cgp::mat9 AqqQuad=cgp::mat9();
    cgp::mat39 ApqQuad;
    cgp::buffer<cgp::vec9> qQuad = {};
    


    

    void initialize(float c = 0.3f, cgp::vec3 globalPosition = cgp::vec3(0.0, 0.0, 0.0),cgp::vec3 anglesEuler=cgp::vec3(0.0f,3.14/4,0.0));  // fill the shape with particles
};
