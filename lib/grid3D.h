#include <vector>

class Grid3D
{
public:
    int Nx, Ny, Nz;

    Grid3D(int nx, int ny, int nz)
        : Nx(nx), Ny(ny), Nz(nz),
          Ex(nx * ny * nz, 0.0),
          Ey(nx * ny * nz, 0.0),
          Ez(nx * ny * nz, 0.0),
          Hx(nx * ny * nz, 0.0),
          Hy(nx * ny * nz, 0.0),
          Hz(nx * ny * nz, 0.0)
    {}

    inline std::size_t id(int i, int j, int k) const
    {
        return static_cast<std::size_t>(
            (k * Ny + j) * Nx + i
        );
    }

    double &ex(int i, int j, int k) { return Ex[id(i,j,k)]; }
    double &ey(int i, int j, int k) { return Ey[id(i,j,k)]; }
    double &ez(int i, int j, int k) { return Ez[id(i,j,k)]; }

    double &hx(int i, int j, int k) { return Hx[id(i,j,k)]; }
    double &hy(int i, int j, int k) { return Hy[id(i,j,k)]; }
    double &hz(int i, int j, int k) { return Hz[id(i,j,k)]; }

private:
    std::vector<double> Ex, Ey, Ez;
    std::vector<double> Hx, Hy, Hz;
};