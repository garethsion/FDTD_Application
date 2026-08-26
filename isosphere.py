import numpy as np
import glob
import plotly.graph_objects as go

# Note: this assumes you have saved the entire 3D Ez field.
# The CSV output from the previous program only contains a 2D slice.

# Example synthetic demonstration:
Nx, Ny, Nz = 60, 60, 60

x = np.linspace(-1, 1, Nx)
y = np.linspace(-1, 1, Ny)
z = np.linspace(-1, 1, Nz)

X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

# Example spherical wave
R = np.sqrt(X**2 + Y**2 + Z**2)
field = np.exp(-20 * (R - 0.5)**2) * np.cos(40 * (R - 0.5))

fig = go.Figure(
    data=go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=field.flatten(),
        isomin=0.3,
        isomax=1.0,
        surface_count=5,
        colorscale="RdBu",
        caps=dict(x_show=False, y_show=False, z_show=False)
    )
)

fig.update_layout(
    title="3D FDTD Electric Field",
    scene=dict(
        xaxis_title="x",
        yaxis_title="y",
        zaxis_title="z"
    )
)

fig.show()