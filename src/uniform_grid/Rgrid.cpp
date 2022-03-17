#include "Rgrid.hpp"

using namespace cgp;

void Rgrid::initialize(int resolution)
{
    m_resolution = resolution;
    grid.clear();
    grid.resize(cgp::int3(resolution, resolution, resolution));
    m_minBox = cgp::vec3(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity());
    m_maxBox = cgp::vec3(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
}

void Rgrid::updateMinMax(cgp::vec3 pos) {
    m_maxBox.x = std::max(m_maxBox.x, pos.x + 0.005f);
    m_minBox.x = std::min(m_minBox.x, pos.x - 0.005f);
    m_maxBox.y = std::max(m_maxBox.y, pos.y + 0.005f);
    m_minBox.y = std::min(m_minBox.y, pos.y - 0.005f);
    m_maxBox.z = std::max(m_maxBox.z, pos.z + 0.005f);
    m_minBox.z = std::min(m_minBox.z, pos.z - 0.005f);
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
