#include "cloth.hpp"

using namespace cgp;


void cloth_structure::initialize(int N_samples_edge_arg)
{
    assert_cgp(N_samples_edge_arg > 3, "N_samples_edge=" + str(N_samples_edge_arg) + " should be > 3");

    position.clear();
    normal.clear();
    velocity.clear();
    force.clear();

    position.resize(N_samples_edge_arg, N_samples_edge_arg);
    normal.resize(N_samples_edge_arg, N_samples_edge_arg);
    velocity.resize(N_samples_edge_arg, N_samples_edge_arg);
    force.resize(N_samples_edge_arg, N_samples_edge_arg);
    

    float const z0 = 1.0f;
    mesh const cloth_mesh = mesh_primitive_grid({ -0.5f,0,z0 }, { 0.5f,0,z0 }, { 0.5f,1,z0 }, { -0.5f,1,z0 }, N_samples_edge_arg, N_samples_edge_arg).fill_empty_field();
    position = grid_2D<vec3>::from_buffer(cloth_mesh.position, N_samples_edge_arg, N_samples_edge_arg);
    normal = grid_2D<vec3>::from_buffer(cloth_mesh.normal, N_samples_edge_arg, N_samples_edge_arg);
    triangle_connectivity = cloth_mesh.connectivity;
}

void cloth_structure::update_normal()
{
    normal_per_vertex(position.data, triangle_connectivity, normal.data);
}

int cloth_structure::N_samples_edge() const
{
    return position.dimension.x;
}



void cloth_structure_drawable::initialize(int N_samples_edge)
{
    mesh const cloth_mesh = mesh_primitive_grid({ 0,0,0 }, {1,0,0 }, { 1,1,0 }, { 0,1,0 }, N_samples_edge, N_samples_edge).fill_empty_field();

    drawable.clear();
    drawable.initialize(cloth_mesh, "Cloth");
    drawable.shading.phong.specular = 0.0f;
}


void cloth_structure_drawable::update(cloth_structure const& cloth)
{    
    drawable.update_position(cloth.position.data);
    drawable.update_normal(cloth.normal.data);
}
