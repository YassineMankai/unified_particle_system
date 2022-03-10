#include "Rgrid.hpp"

using namespace cgp;

void Rgrid::initialize(cgp::vec3 minBox, cgp::vec3 maxBox, int resolution)
{
    m_resolution = resolution;
    grid.clear();
    grid.resize(cgp::int3(resolution, resolution, resolution));
    m_minBox = minBox;
    m_maxBox = maxBox;
}

cgp::int3 Rgrid::posToIndex(cgp::vec3 pos) const
{
    cgp::vec3 posNormalized = (m_resolution - 1) * (pos - m_minBox) / (m_maxBox - m_minBox);
    return cgp::int3(static_cast<int>(posNormalized.x), static_cast<int>(posNormalized.y), static_cast<int>(posNormalized.z));
}

void Rgrid::insert(vec3 pos, int index)
{
    cgp::int3 cell = posToIndex(pos);
    grid[cell].push_back(index);
}

cgp::buffer<int> Rgrid::getNeighborhood(vec3 pos) const
{
    cgp::int3 cell = posToIndex(pos);
    cgp::buffer<int> res;

    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                cgp::int3 new_cell = cell + cgp::int3(dx, dy, dz);
                if (new_cell.x >= 0 && new_cell.x < m_resolution
                 && new_cell.y >= 0 && new_cell.y < m_resolution
                 && new_cell.z >= 0 && new_cell.z < m_resolution)
                {
                    res.push_back(grid[new_cell]);
                }
            }
        }
    }
    return res;
}
