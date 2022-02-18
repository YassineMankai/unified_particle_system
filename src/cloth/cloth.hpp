#pragma once

#include "cgp/cgp.hpp"

// Stores the buffers representing the cloth vertices
struct cloth_structure
{    
    // Buffers are stored as 2D grid that can be accessed as grid(ku,kv)
    cgp::grid_2D<cgp::vec3> position;  
    cgp::grid_2D<cgp::vec3> velocity;  
    cgp::grid_2D<cgp::vec3> force;
    cgp::grid_2D<cgp::vec3> normal;

    // Also stores the triangle connectivity used to update the normals
    cgp::buffer<cgp::uint3> triangle_connectivity;

    
    void initialize(int N_samples_edge);  // Initialize a square flat cloth
    void update_normal();       // Call this function every time the cloth is updated before its draw
    int N_samples_edge() const; // Number of vertex along one dimension of the grid
};


// Helper structure and functions to draw a cloth
// ********************************************** //
struct cloth_structure_drawable
{
    cgp::mesh_drawable drawable;

    void initialize(int N_sample_edge);
    void update(cloth_structure const& cloth);
};

template <typename ENVIRONMENT>
void draw(cloth_structure_drawable const& drawable, ENVIRONMENT const& environment)
{
    draw(drawable.drawable, environment);
}
template <typename ENVIRONMENT>
void draw_wireframe(cloth_structure_drawable const& drawable, ENVIRONMENT const& environment)
{
    draw_wireframe(drawable.drawable, environment);
}