#pragma once

#include "cgp/cgp.hpp"
#include "cgp/containers/grid/grid_3D/grid_3D.hpp"
using namespace cgp;


class Rgrid
{
private:
    int m_resolution;
    cgp::grid_3D<buffer<int>> grid;
    cgp::vec3 m_minBox;
    cgp::vec3 m_maxBox;
    cgp::int3 posToIndex(cgp::vec3 pos) const;

public:
    void initialize(cgp::vec3 minBox, cgp::vec3 maxBox, int resolution);
    void insert(vec3 pos, int index);
    cgp::buffer<int> getNeighborhood(vec3 pos) const;
};
