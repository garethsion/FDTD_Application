class SimplePML {
    // Simple lossy boundary coefficient.
    // This is NOT a true PML.
    // It damps fields close to the outer boundary to reduce reflections.
    public:
        double dampingCoefficient(int i, int j, int k, int Nx, int Ny, int Nz, int boundaryWidth)
        {
            int d = std::min({i, j, k, Nx - 1 - i, Ny - 1 - j, Nz - 1 - k});

            if (d >= boundaryWidth)
                return 1.0;

            double x = static_cast<double>(boundaryWidth - d)
                    / static_cast<double>(boundaryWidth);

            // Smooth polynomial damping.
            return std::exp(-0.15 * x * x);
        }
};