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
OFFSET_KERNEL_V = np.array([-5.0, 0.0, 0.0])
OFFSET_KERNEL_G = np.array([5.0, 0.0, -vg.CFG.param_stk_dist_mu])

# ------------------------------------------------------------------------------
def stamp_line_v(grid: vg.Grid, angle: float):
    vdir = np.array([np.cos(angle), np.sin(angle), 0.0])
    kernel = vg.KernelCylinder(
        radius = RADIUS_KERNEL_V,
        vdirection = vdir,
        width = WIDTH_KERNEL_V,
        deltas = grid.box.deltas,
        dtype = float,
    )
    center_stamp_at = RADIUS_KERNEL_V * vdir + OFFSET_KERNEL_V
    kernel.stamp(grid, center_stamp_at)


# ------------------------------------------------------------------------------
def main():
    ##### PART 0: INIT GRID
    box = vg.Box(
        origin = np.array([-20.0, -20.0, -20.0]),
        resolution = np.array([160, 160, 160]),
        deltas = np.array([0.25, 0.25, 0.25]),
    )
    grid = vg.Grid(box)
    print(grid.arr.shape, box.min_coords, box.max_coords)

    ##### PART 1: DRAW THE V
    angle_v = (3/8)*np.pi # must be less than pi/2
    stamp_line_v(grid, angle_v)
    stamp_line_v(grid, np.pi - angle_v)


    ##### PART 2: DRAW THE G


    ##### PART 3: SAVE
    grid.save(PATH_MRC_OUT)


################################################################################
if __name__ == "__main__":
    PATH_MRC_OUT = Path("examples/logo/logo.mrc")
    main()


################################################################################
# python3 examples/logo/gen_logo_mrc.py
