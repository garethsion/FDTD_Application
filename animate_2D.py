import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

files = sorted(
    glob.glob("ez_*.csv"),
    key=lambda f: int(f.split("_")[1].split(".")[0])
)

if not files:
    raise RuntimeError("No ez_*.csv files found.")

# Read first frame to determine grid size
first = np.loadtxt(files[0], delimiter=",")
ny, nx = first.shape

fig, ax = plt.subplots(figsize=(8, 7))

image = ax.imshow(
    first,
    origin="lower",
    extent=[0, nx, 0, ny],
    interpolation="bilinear",
    cmap="RdBu",
    vmin=-1,
    vmax=1
)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("3D FDTD: Ez field")

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("Ez")

def update(frame):
    data = np.loadtxt(files[frame], delimiter=",")

    image.set_data(data)

    step = files[frame].split("_")[1].split(".")[0]
    ax.set_title(f"3D FDTD: Ez field — step {step}")

    return [image]

animation = FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=50,
    blit=True
)

plt.show()