import numpy as np

### simulate having "volgrids" installed as a package
### this way it's not necessary to install the repo to run this script
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
### you can remove the previous two lines if volgrids is installed

import volgrids as vg

WIDTH_KERNEL_V = 1.0
RADIUS_KERNEL_V = 5.0
RADIUS_KERNEL_G = 5.0
OFFSET_KERNEL_V = np.array([-2.0, 0.0, 0.0])
OFFSET_KERNEL_G = np.array([3.0, -1.0, -RADIUS_KERNEL_G])

# ------------------------------------------------------------------------------
def stamp_line_v(grid: vg.Grid, angle: float):
    vdir = np.array([np.cos(angle), np.sin(angle), 0.0])
    kernel = vg.KernelCylinder(
        radius = RADIUS_KERNEL_V,
        vnormal = vdir,
        width_inner = 0, width_outer = WIDTH_KERNEL_V,
        height = RADIUS_KERNEL_V,
        deltas = grid.box.deltas,
    )
    center_stamp_at = RADIUS_KERNEL_V * vdir + OFFSET_KERNEL_V
    kernel.stamp(grid, center_stamp_at)


# ------------------------------------------------------------------------------
def main():
    ##### PART 0: INIT GRID
    box = vg.Box(
        origin = np.array([-10.0, -10.0, -10.0]),
        resolution = np.array([250, 250, 250]),
        deltas = np.array([0.10, 0.10, 0.10]),
    )
    grid = vg.Grid(box)
    print(grid.arr.shape, box.min_coords, box.max_coords)


    ##### PART 1: DRAW THE V
    angle_v = (3/8)*np.pi # must be less than pi/2
    stamp_line_v(grid, angle_v)
    stamp_line_v(grid, np.pi - angle_v)


    ##### PART 2: DRAW THE G
    normal = np.array([0.0, 0.0, 1.0])
    width_inner_g = RADIUS_KERNEL_G - 2*WIDTH_KERNEL_V
    width_outer_g = RADIUS_KERNEL_G

    kernel = vg.KernelCylinder(
        radius = RADIUS_KERNEL_G,
        deltas = box.deltas,
        vnormal = normal,
        width_inner = width_inner_g,
        width_outer = width_outer_g,
        height = WIDTH_KERNEL_V,
    )
    center_stamp_at = RADIUS_KERNEL_G + OFFSET_KERNEL_G
    kernel.stamp(grid, center_stamp_at)


    kernel = vg.Kernel(
        RADIUS_KERNEL_G, box.deltas,
        dtype = bool, kop = vg.KOperation.MASK
    )
    center_stamp_at = RADIUS_KERNEL_G + OFFSET_KERNEL_G
    center_stamp_at[:2] += RADIUS_KERNEL_G
    kernel.stamp(grid, center_stamp_at)


    kernel = vg.KernelCylinder(
        radius = RADIUS_KERNEL_G,
        deltas = box.deltas,
        vnormal = np.array([1.0, 0.0, 0.0]),
        width_inner = 0,
        width_outer = WIDTH_KERNEL_V,
        height = width_outer_g/2*0.90, # small adjustment
    )
    center_stamp_at = RADIUS_KERNEL_G + OFFSET_KERNEL_G
    center_stamp_at[0] += RADIUS_KERNEL_G/2
    center_stamp_at[1] -= RADIUS_KERNEL_G*0.05 # small adjustment
    kernel.stamp(grid, center_stamp_at)

    ##### PART 3: SAVE
    grid.save(PATH_MRC_OUT)


################################################################################
if __name__ == "__main__":
    PATH_MRC_OUT = Path("examples/logo/logo.mrc")
    main()


################################################################################
# python3 examples/logo/gen_logo_mrc.py
