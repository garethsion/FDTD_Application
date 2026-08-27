#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include "sources.h"
#include "simple_pml.h"
#include "grid3D.h"
#include "materials.h"

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0e-7 * M_PI;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

int main()
{
    // Simulation parameters
    const int Nx = 120;
    const int Ny = 120;
    const int Nz = 120;

    const double dx = 1.0e-3;
    const double dy = dx;
    const double dz = dx;

    const double dt = 0.95 / (c0 * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
    const int nSteps = 1000;

    std::cout << "dt = " << dt << " s\n";
    std::cout << "Steps = " << nSteps << "\n";

    // Grid
    Grid3D grid(Nx, Ny, Nz);

    // Material
    IdealMaterial material(Nx, Ny, Nz, 4.0, 1.0);

    // Source
    const int sx = Nx / 2;
    const int sy = Ny / 2;
    const int sz = Nz / 2;

    const double tau = 20.0 * dt;
    const double t0 = 6.0 * tau;

    // Boundary damping region
    const int boundaryWidth = 12;

    // Save a 2D Ez slice every N steps
    const int outputEvery = 20;

    // FDTD time stepping
    for (int n = 0; n < nSteps; ++n) {
        const double t = n * dt;

        // =====================================================================
        // Update H fields
        // H(n+1/2) = H(n-1/2) - dt/mu * curl(E)
        // =====================================================================

        for (int k = 0; k < Nz - 1; ++k) {
            for (int j = 0; j < Ny - 1; ++j) {
                for (int i = 0; i < Nx - 1; ++i) {
                    const std::size_t p = material.matId(i, j, k);
                    const double mu = mu0 * material.muR[p];

                    // curl(E)_x
                    const double dEz_dy = (grid.ez(i, j + 1, k) - grid.ez(i, j, k)) / dy;
                    const double dEy_dz = (grid.ey(i, j, k + 1) - grid.ey(i, j, k)) / dz;
                    grid.hx(i, j, k) -= (dt / mu) * (dEz_dy - dEy_dz);

                    // curl(E)_y
                    const double dEx_dz = (grid.ex(i, j, k + 1) - grid.ex(i, j, k)) / dz;
                    const double dEz_dx = (grid.ez(i + 1, j, k) - grid.ez(i, j, k)) / dx;
                    grid.hy(i, j, k) -= (dt / mu) * (dEx_dz - dEz_dx);

                    // curl(E)_z
                    const double dEy_dx = (grid.ey(i + 1, j, k) - grid.ey(i, j, k)) / dx;
                    const double dEx_dy = (grid.ex(i, j + 1, k) - grid.ex(i, j, k)) / dy;
                    grid.hz(i, j, k) -= (dt / mu) * (dEy_dx - dEx_dy);
                }
            }
        }

        // =====================================================================
        // Update E fields
        // E(n+1) = E(n) + dt/eps * curl(H)
        // =====================================================================

        for (int k = 1; k < Nz - 1; ++k) {
            for (int j = 1; j < Ny - 1; ++j) {
                for (int i = 1; i < Nx - 1; ++i) {
                    const std::size_t p = material.matId(i, j, k);
                    const double eps = eps0 * material.epsR[p];

                    // curl(H)_x
                    const double dHz_dy = (grid.hz(i, j, k) - grid.hz(i, j - 1, k)) / dy;
                    const double dHy_dz = (grid.hy(i, j, k) - grid.hy(i, j, k - 1)) / dz;
                    grid.ex(i, j, k) += (dt / eps) * (dHz_dy - dHy_dz);

                    // curl(H)_y
                    const double dHx_dz = (grid.hx(i, j, k) - grid.hx(i, j, k - 1)) / dz;
                    const double dHz_dx = (grid.hz(i, j, k) - grid.hz(i - 1, j, k)) / dx;
                    grid.ey(i, j, k) += (dt / eps) * (dHx_dz - dHz_dx);

                    // curl(H)_z
                    const double dHy_dx = (grid.hy(i, j, k) - grid.hy(i - 1, j, k)) / dx;
                    const double dHx_dy = (grid.hx(i, j, k) - grid.hx(i, j - 1, k)) / dy;
                    grid.ez(i, j, k) += (dt / eps) * (dHy_dx - dHx_dy);
                }
            }
        }

        // Inject source //
        grid.ez(sx, sy, sz) += Sources().gaussianPulse(t, t0, tau);

        // Apply PML //
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    const double d = SimplePML().dampingCoefficient(i, j, k, Nx, Ny, Nz, boundaryWidth);

                    grid.ex(i, j, k) *= d;
                    grid.ey(i, j, k) *= d;
                    grid.ez(i, j, k) *= d;

                    grid.hx(i, j, k) *= d;
                    grid.hy(i, j, k) *= d;
                    grid.hz(i, j, k) *= d;
                }
            }
        }

        // =====================================================================
        // Save central z-slice
        // =====================================================================

        if (n % outputEvery == 0) {
            const std::string filename = "outputs/ez_" + std::to_string(n) + ".csv";
            std::ofstream file(filename);
            const int kz = Nz / 2;

            if (!file) {
                std::cerr << "Could not open output file: " << filename << "\n";
                return 1;
            }

            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    file << grid.ez(i, j, kz);

                    if (i + 1 < Nx) {
                        file << ',';
                    }
                }

                file << '\n';
            }

            file.close();
            std::cout << "Step " << n << " / " << nSteps << "\n";
        }
    }
    std::cout << "Simulation complete.\n";
    return 0;
}