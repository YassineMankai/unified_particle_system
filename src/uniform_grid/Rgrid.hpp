#pragma once

#include "cgp/cgp.hpp"
#include "cgp/containers/grid/grid_3D/grid_3D.hpp"
using namespace cgp;


struct Rgrid
{
private:
    int m_resolution;
    cgp::grid_3D<buffer<int>> grid;
    cgp::vec3 m_minBox;
    cgp::vec3 m_maxBox;
    cgp::int3 posToIndex(cgp::vec3 pos) const;

public:
    void initialize(int resolution);
    void insert(vec3 pos, int index);
    std::vector<int> getNeighborhood(vec3 pos) const;
    void updateMinMax(cgp::vec3 pos);
};
