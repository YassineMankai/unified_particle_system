#pragma once

#include "custom_special_types.hpp"
#include "cgp/cgp.hpp"



// SPH Particle
struct particle_element
{
    cgp::vec3 position; // Position
    cgp::vec3 velocity; // Speed
    cgp::vec3 force; // Force
    int nbConstraint = 0;
    int nbConstraintTest = 0;
    cgp::vec3 dx = cgp::vec3(0,0,0); 
    cgp::vec3 dx_friction_and_restitution = cgp::vec3(0,0,0); 
    float mass = 1;
    int phase;

    particle_element() : position{ 0,0,0 }, velocity{ 0,0,0 }, force{ 0,0,0 } {}
};

// Scene creation parameters
struct scene_parameters
{
    float shape_size = 0.06f;
};

enum ShapeType {QUADRATIC, RIGID, CLOTH};

// Stores the buffers representing the shape vertices
struct shape_structure
{    
    ShapeType type;
    cgp::vec3 com0; //intial center of mass
    int relativeLocationsOffset;
    int nbParticles;
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
    
    //cloth variable
    int height = 0;
    int width = 0;
    float bendLength0;
    float structLength0;
    float shearLength0;
    int getGlobalIndex(int i, int j);
    int getNeighbour(int p_Index, int i,int j);
    

};
