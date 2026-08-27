// // ============================================================
// // CPML coefficients for one spatial coordinate
// //
// // We need coefficients at:
// //
// //   integer positions : electric-field update locations
// //   half positions    : magnetic-field update locations
// //
// // Each set contains:
// //
// //   kappa
// //   b
// //   c
// //
// // The CPML modified derivative is:
// //
// //   d~f/dx = (1/kappa) df/dx + psi
// //
// // with:
// //
// //   psi(n+1) = b psi(n) + c df/dx
// // ============================================================

// struct CPML1D
// {
//     std::vector<double> kappaInt;
//     std::vector<double> bInt;
//     std::vector<double> cInt;

//     std::vector<double> kappaHalf;
//     std::vector<double> bHalf;
//     std::vector<double> cHalf;
// };


// double sigmaProfile(
//     double distance,
//     double thickness,
//     double sigmaMax,
//     int order)
// {
//     if (distance >= thickness)
//         return 0.0;

//     const double u =
//         (thickness - distance) / thickness;

//     return sigmaMax * std::pow(u, order);
// }


// double kappaProfile(
//     double distance,
//     double thickness,
//     double kappaMax,
//     int order)
// {
//     if (distance >= thickness)
//         return 1.0;

//     const double u =
//         (thickness - distance) / thickness;

//     return 1.0 +
//            (kappaMax - 1.0) *
//            std::pow(u, order);
// }


// double alphaProfile(
//     double distance,
//     double thickness,
//     double alphaMax)
// {
//     if (distance >= thickness)
//         return 0.0;

//     const double u =
//         (thickness - distance) / thickness;

//     // Maximum at the inner edge of the PML,
//     // decreasing toward the outer boundary.
//     return alphaMax * (1.0 - u);
// }


// void makeCPMLCoefficients(
//     int N,
//     double dx,
//     int nPML,
//     double targetReflection,
//     int polynomialOrder,
//     double kappaMax,
//     double alphaMax,
//     CPML1D& cpml)
// {
//     cpml.kappaInt.resize(N + 1);
//     cpml.bInt.resize(N + 1);
//     cpml.cInt.resize(N + 1);

//     cpml.kappaHalf.resize(N);
//     cpml.bHalf.resize(N);
//     cpml.cHalf.resize(N);

//     const double thickness =
//         static_cast<double>(nPML) * dx;

//     // Standard PML conductivity scaling.
//     const double sigmaMax =
//         -(polynomialOrder + 1.0) *
//         std::log(targetReflection) /
//         (2.0 * eta0 * thickness);

//     // --------------------------------------------------------
//     // Integer positions.
//     //
//     // Coordinates:
//     //
//     //   x = i * dx
//     //
//     // Boundaries are at 0 and N*dx.
//     // --------------------------------------------------------

//     for (int i = 0; i <= N; ++i)
//     {
//         const double distance =
//             std::min(
//                 static_cast<double>(i) * dx,
//                 static_cast<double>(N - i) * dx
//             );

//         const double sigma =
//             sigmaProfile(
//                 distance,
//                 thickness,
//                 sigmaMax,
//                 polynomialOrder
//             );

//         const double kappa =
//             kappaProfile(
//                 distance,
//                 thickness,
//                 kappaMax,
//                 polynomialOrder
//             );

//         const double alpha =
//             alphaProfile(
//                 distance,
//                 thickness,
//                 alphaMax
//             );

//         // This is filled later by setTimeStep().
//         cpml.kappaInt[i] = kappa;

//         // Store temporary b/c information using
//         // negative sentinel values is ugly, so initially
//         // use sigma/alpha encoded through reconstruction.
//         //
//         // Instead, calculate with a temporary local formula
//         // in setTimeStep().
//         cpml.bInt[i] = sigma;
//         cpml.cInt[i] = alpha;
//     }

//     // --------------------------------------------------------
//     // Half positions.
//     //
//     // Coordinates:
//     //
//     //   x = (i + 0.5) * dx
//     //
//     // --------------------------------------------------------

//     for (int i = 0; i < N; ++i)
//     {
//         const double x =
//             (static_cast<double>(i) + 0.5) * dx;

//         const double distance =
//             std::min(
//                 x,
//                 static_cast<double>(N) * dx - x
//             );

//         const double sigma =
//             sigmaProfile(
//                 distance,
//                 thickness,
//                 sigmaMax,
//                 polynomialOrder
//             );

//         const double kappa =
//             kappaProfile(
//                 distance,
//                 thickness,
//                 kappaMax,
//                 polynomialOrder
//             );

//         const double alpha =
//             alphaProfile(
//                 distance,
//                 thickness,
//                 alphaMax
//             );

//         cpml.kappaHalf[i] = kappa;
//         cpml.bHalf[i] = sigma;
//         cpml.cHalf[i] = alpha;
//     }
// }


// // ============================================================
// // Finalize CPML coefficients after dt is known.
// //
// // Input bInt/cInt temporarily contain:
// //
// //   bInt = sigma
// //   cInt = alpha
// //
// // and similarly for half positions.
// // ============================================================

// void finalizeCPML(
//     CPML1D& cpml,
//     double dt)
// {
//     auto convert = [dt](double kappa,
//                         double sigma,
//                         double alpha)
//     {
//         // sigma/epsilon has units of 1/s.
//         const double s =
//             sigma / eps0;

//         const double b =
//             std::exp(
//                 -(s / kappa + alpha) * dt
//             );

//         double c = 0.0;

//         const double denominator =
//             s * kappa +
//             alpha * kappa * kappa;

//         if (std::abs(denominator) > 1e-30)
//         {
//             c =
//                 s / denominator *
//                 (b - 1.0);
//         }

//         return std::pair<double, double>{b, c};
//     };

//     for (std::size_t i = 0;
//          i < cpml.kappaInt.size();
//          ++i)
//     {
//         const double kappa =
//             cpml.kappaInt[i];

//         const double sigma =
//             cpml.bInt[i];

//         const double alpha =
//             cpml.cInt[i];

//         auto [b, c] =
//             convert(kappa, sigma, alpha);

//         cpml.bInt[i] = b;
//         cpml.cInt[i] = c;
//     }

//     for (std::size_t i = 0;
//          i < cpml.kappaHalf.size();
//          ++i)
//     {
//         const double kappa =
//             cpml.kappaHalf[i];

//         const double sigma =
//             cpml.bHalf[i];

//         const double alpha =
//             cpml.cHalf[i];

//         auto [b, c] =
//             convert(kappa, sigma, alpha);

//         cpml.bHalf[i] = b;
//         cpml.cHalf[i] = c;
//     }
// }