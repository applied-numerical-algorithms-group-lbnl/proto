#include <cmath>
#include <vector>
#include <array>

int linearIndex(int i, int j, int k, int N)
{
    return N*N*k + N*j + i;
}

std::array<int, 3> getLowCorner(float x, float y, float z, int N)
{
    int i0 = static_cast<int>(std::floor(x));
    int j0 = static_cast<int>(std::floor(y));
    int k0 = static_cast<int>(std::floor(z));
    i0 = std::max(0, std::min(i0, N - 2));
    j0 = std::max(0, std::min(j0, N - 2));
    k0 = std::max(0, std::min(k0, N - 2));
    return {i0, j0, k0};
}

float gridToImplicit(float x, float y, float z, const std::vector<float>& sdfGrid, int N) {
    // Convert world coordinates to grid coordinates
    // Assuming inputs are in unit cube
    float gridX = x * (N - 1);
    float gridY = y * (N - 1);
    float gridZ = z * (N - 1);
    
    // get low corner indices
    int i0 = static_cast<int>(std::floor(gridX));
    int j0 = static_cast<int>(std::floor(gridY));
    int k0 = static_cast<int>(std::floor(gridZ));
    
    // bound check
    i0 = std::max(0, std::min(i0, N - 2));
    j0 = std::max(0, std::min(j0, N - 2));
    k0 = std::max(0, std::min(k0, N - 2));
    
    auto 

    // get upper corner indices
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;
    
    // Calculate interpolation weights
    float tx = gridX - i0;
    float ty = gridY - j0;
    float tz = gridZ - k0;
    
    // Fetch the 8 values from the grid
    float v000 = sdfGrid[linearIndex(i0, j0, k0, N)];
    float v001 = sdfGrid[linearIndex(i0, j0, k1, N)];
    float v010 = sdfGrid[linearIndex(i0, j1, k0, N)];
    float v011 = sdfGrid[linearIndex(i0, j1, k1, N)];
    float v100 = sdfGrid[linearIndex(i1, j0, k0, N)];
    float v101 = sdfGrid[linearIndex(i1, j0, k1, N)];
    float v110 = sdfGrid[linearIndex(i1, j1, k0, N)];
    float v111 = sdfGrid[linearIndex(i1, j1, k1, N)];
    
    // Perform trilinear interpolation
    float v00 = v000 * (1 - tx) + v100 * tx;
    float v01 = v001 * (1 - tx) + v101 * tx;
    float v10 = v010 * (1 - tx) + v110 * tx;
    float v11 = v011 * (1 - tx) + v111 * tx;
    
    float v0 = v00 * (1 - ty) + v10 * ty;
    float v1 = v01 * (1 - ty) + v11 * ty;
    
    float v = v0 * (1 - tz) + v1 * tz;
    
    return v;
}