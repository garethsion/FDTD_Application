#include <cstddef>
#include <vector>

class IdealMaterial
{
public:
    int Nx, Ny, Nz;

    std::vector<double> epsR;
    std::vector<double> muR;

    IdealMaterial(int nx, int ny, int nz, double eps_r = 1.0, double mu_r = 1.0)
        : Nx(nx),
          Ny(ny),
          Nz(nz),
          epsR(static_cast<std::size_t>(nx) * ny * nz, 1.0),
          muR(static_cast<std::size_t>(nx) * ny * nz, 1.0)
    {
        for (int k = 45; k < 75; ++k) {
            for (int j = 45; j < 75; ++j) {
                for (int i = 45; i < 75; ++i) {
                    epsR[matId(i, j, k)] = eps_r;
                    muR[matId(i, j, k)] = mu_r;
                }
            }
        }
    }

    std::size_t matId(int i, int j, int k) const
    {
        return static_cast<std::size_t>((k * Ny + j) * Nx + i);
    }
};