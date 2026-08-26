#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

class Sources
{
public:
    /**
     * @brief Evaluate a Gaussian pulse.
     *
     * g(t) = exp(-((t - t0) / tau)^2)
     *
     * @param t   Current time [s].
     * @param t0  Time at which the pulse reaches its maximum [s].
     * @param tau Pulse-width parameter [s].
     *
     * @return Gaussian pulse amplitude.
     */
    static double gaussianPulse(double t, double t0, double tau)
    {
        if (tau <= 0.0) {
            throw std::invalid_argument(
                "gaussianPulse: tau must be greater than zero"
            );
        }

        const double x = (t - t0) / tau;
        return std::exp(-x * x);
    }


    /**
     * @brief Electric dipole source.
     *
     * Represents
     *
     *      p(t) = p0 * direction * g(t)
     *
     * where p0 is the peak dipole moment and g(t) is a Gaussian pulse.
     *
     * @param t         Current time [s].
     * @param t0        Centre time of the Gaussian pulse [s].
     * @param tau       Pulse-width parameter [s].
     * @param p0        Peak electric dipole moment [C m].
     * @param direction Dipole orientation vector.
     *
     * @return Dipole moment vector [C m].
     */
    static std::array<double, 3> electricDipole(
        double t,
        double t0,
        double tau,
        double p0,
        const std::array<double, 3>& direction)
    {
        // Normalize dipole direction
        const double norm = std::sqrt(
            direction[0] * direction[0] +
            direction[1] * direction[1] +
            direction[2] * direction[2]
        );

        if (norm == 0.0) {
            throw std::invalid_argument(
                "electricDipole: direction cannot be zero"
            );
        }

        const double pulse = gaussianPulse(t, t0, tau);

        const double amplitude = p0 * pulse / norm;

        return {
            amplitude * direction[0],
            amplitude * direction[1],
            amplitude * direction[2]
        };
    }
};