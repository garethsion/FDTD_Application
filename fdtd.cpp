#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

constexpr double c0   = 299792458.0;
constexpr double mu0  = 4.0e-7 * M_PI;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

// -----------------------------------------------------------------------------
// 3D FDTD simulation using the Yee algorithm.
//
// Grid convention:
//   Ex(i,j,k)
//   Ey(i,j,k)
//   Ez(i,j,k)
//   Hx(i,j,k)
//   Hy(i,j,k)
//   Hz(i,j,k)
//
// For simplicity, all components are stored on the same integer grid.
// This is slightly less memory-efficient than a staggered Yee implementation,
// but the update equations below are arranged to approximate the Yee scheme.
//
// The code injects a Gaussian pulse into Ez and writes Ez slices to CSV.
// -----------------------------------------------------------------------------

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


// Gaussian pulse
double gaussianPulse(double t, double t0, double tau)
{
    const double x = (t - t0) / tau;
    return std::exp(-x * x);
}


// Simple lossy boundary coefficient.
//
// This is NOT a true PML.
// It damps fields close to the outer boundary to reduce reflections.
double dampingCoefficient(
    int i, int j, int k,
    int Nx, int Ny, int Nz,
    int boundaryWidth)
{
    int d = std::min({
        i,
        j,
        k,
        Nx - 1 - i,
        Ny - 1 - j,
        Nz - 1 - k
    });

    if (d >= boundaryWidth)
        return 1.0;

    double x = static_cast<double>(boundaryWidth - d)
             / static_cast<double>(boundaryWidth);

    // Smooth polynomial damping.
    return std::exp(-0.15 * x * x);
}


int main()
{
    // -------------------------------------------------------------------------
    // Simulation parameters
    // -------------------------------------------------------------------------

    const int Nx = 120;
    const int Ny = 120;
    const int Nz = 120;

    const double dx = 1.0e-3;
    const double dy = dx;
    const double dz = dx;

    // CFL-stable timestep in 3D:
    //
    // dt <= 1 / [ c0 * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2) ]
    //
    const double dt =
        0.95 /
        (c0 * std::sqrt(
            1.0/(dx*dx) +
            1.0/(dy*dy) +
            1.0/(dz*dz)
        ));

    const int nSteps = 1000;

    std::cout << "dt = " << dt << " s\n";
    std::cout << "Steps = " << nSteps << "\n";

    Grid3D grid(Nx, Ny, Nz);

    // -------------------------------------------------------------------------
    // Material properties
    //
    // Here everything is vacuum.
    // Change epsR/muR to model materials.
    // -------------------------------------------------------------------------

    std::vector<double> epsR(
        static_cast<std::size_t>(Nx) * Ny * Nz, 1.0
    );

    std::vector<double> muR(
        static_cast<std::size_t>(Nx) * Ny * Nz, 1.0
    );

    auto matId = [=](int i, int j, int k)
    {
        return static_cast<std::size_t>(
            (k * Ny + j) * Nx + i
        );
    };

    // Example dielectric cube.
    //
    // Relative permittivity = 4
    //
    for (int k = 45; k < 75; ++k)
    {
        for (int j = 45; j < 75; ++j)
        {
            for (int i = 45; i < 75; ++i)
            {
                epsR[matId(i,j,k)] = 4.0;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Source parameters
    // -------------------------------------------------------------------------

    const int sx = Nx / 2;
    const int sy = Ny / 2;
    const int sz = Nz / 2;

    const double tau = 20.0 * dt;
    const double t0  = 6.0 * tau;

    // Boundary damping region.
    const int boundaryWidth = 12;

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------

    // Save a 2D Ez slice every N steps.
    const int outputEvery = 20;

    // -------------------------------------------------------------------------
    // FDTD time stepping
    // -------------------------------------------------------------------------

    for (int n = 0; n < nSteps; ++n)
    {
        const double t = n * dt;

        // =====================================================================
        // Update H fields
        //
        // Hx(n+1/2) = Hx(n-1/2)
        //            - dt/mu * curl(E)
        //
        // =====================================================================

        for (int k = 0; k < Nz - 1; ++k)
        {
            for (int j = 0; j < Ny - 1; ++j)
            {
                for (int i = 0; i < Nx - 1; ++i)
                {
                    const std::size_t p = matId(i,j,k);

                    const double mu = mu0 * muR[p];

                    // curl(E)_x
                    const double dEz_dy =
                        (grid.ez(i, j + 1, k) -
                         grid.ez(i, j,     k)) / dy;

                    const double dEy_dz =
                        (grid.ey(i, j, k + 1) -
                         grid.ey(i, j, k)) / dz;

                    grid.hx(i,j,k) -=
                        (dt / mu) * (dEz_dy - dEy_dz);

                    // curl(E)_y
                    const double dEx_dz =
                        (grid.ex(i, j, k + 1) -
                         grid.ex(i, j, k)) / dz;

                    const double dEz_dx =
                        (grid.ez(i + 1, j, k) -
                         grid.ez(i, j, k)) / dx;

                    grid.hy(i,j,k) -=
                        (dt / mu) * (dEx_dz - dEz_dx);

                    // curl(E)_z
                    const double dEy_dx =
                        (grid.ey(i + 1, j, k) -
                         grid.ey(i, j, k)) / dx;

                    const double dEx_dy =
                        (grid.ex(i, j + 1, k) -
                         grid.ex(i, j, k)) / dy;

                    grid.hz(i,j,k) -=
                        (dt / mu) * (dEy_dx - dEx_dy);
                }
            }
        }

        // =====================================================================
        // Update E fields
        //
        // E(n+1) = E(n) + dt/eps * curl(H)
        //
        // =====================================================================

        for (int k = 1; k < Nz - 1; ++k)
        {
            for (int j = 1; j < Ny - 1; ++j)
            {
                for (int i = 1; i < Nx - 1; ++i)
                {
                    const std::size_t p = matId(i,j,k);

                    const double eps =
                        eps0 * epsR[p];

                    // curl(H)_x
                    const double dHz_dy =
                        (grid.hz(i, j,     k) -
                         grid.hz(i, j - 1, k)) / dy;

                    const double dHy_dz =
                        (grid.hy(i, j, k) -
                         grid.hy(i, j, k - 1)) / dz;

                    grid.ex(i,j,k) +=
                        (dt / eps) * (dHz_dy - dHy_dz);

                    // curl(H)_y
                    const double dHx_dz =
                        (grid.hx(i, j, k) -
                         grid.hx(i, j, k - 1)) / dz;

                    const double dHz_dx =
                        (grid.hz(i, j, k) -
                         grid.hz(i - 1, j, k)) / dx;

                    grid.ey(i,j,k) +=
                        (dt / eps) * (dHx_dz - dHz_dx);

                    // curl(H)_z
                    const double dHy_dx =
                        (grid.hy(i, j, k) -
                         grid.hy(i - 1, j, k)) / dx;

                    const double dHx_dy =
                        (grid.hx(i, j, k) -
                         grid.hx(i, j - 1, k)) / dy;

                    grid.ez(i,j,k) +=
                        (dt / eps) * (dHy_dx - dHx_dy);
                }
            }
        }

        // =====================================================================
        // Gaussian source
        //
        // Inject into Ez at the center of the grid.
        // =====================================================================

        grid.ez(sx, sy, sz) +=
            gaussianPulse(t, t0, tau);

        // =====================================================================
        // Simple absorbing boundary
        //
        // Again, this is only an approximate absorber.
        // =====================================================================

        for (int k = 0; k < Nz; ++k)
        {
            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    const double d =
                        dampingCoefficient(
                            i, j, k,
                            Nx, Ny, Nz,
                            boundaryWidth
                        );

                    grid.ex(i,j,k) *= d;
                    grid.ey(i,j,k) *= d;
                    grid.ez(i,j,k) *= d;

                    grid.hx(i,j,k) *= d;
                    grid.hy(i,j,k) *= d;
                    grid.hz(i,j,k) *= d;
                }
            }
        }

        // =====================================================================
        // Save a central z-slice
        // =====================================================================

        if (n % outputEvery == 0)
        {
            const std::string filename =
                "ez_" + std::to_string(n) + ".csv";

            std::ofstream file(filename);

            const int kz = Nz / 2;

            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    file << grid.ez(i,j,kz);

                    if (i + 1 < Nx)
                        file << ',';
                }

                file << '\n';
            }

            file.close();

            std::cout
                << "Step " << n
                << " / " << nSteps
                << "\n";
        }
    }

    std::cout << "Simulation complete.\n";

    return 0;
}