import numpy as np

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelSphere(vg.Kernel):
    """For generating simple boolean spheres (e.g. for masks)"""
    def __init__(self,
        radius, deltas,
        dtype = bool, kop: vg.KOperation = vg.KOperation.OR,
    ):
        super().__init__(radius, deltas, dtype, kop)
        self.arr[self.dists < radius] = 1


# //////////////////////////////////////////////////////////////////////////////
class KernelCylinder(vg.Kernel):
    """For generating boolean cylinders, disks and rings"""
    def __init__(self,
        radius, deltas, vnormal: np.ndarray,
        width_inner: float, width_outer: float, height: float,
        dtype = bool, kop: vg.KOperation = vg.KOperation.OR,
    ):
        super().__init__(radius, deltas, dtype, kop)
        self.arr.fill(1)

        w = vg.Math.get_projection_height(self.coords, vnormal)
        mask_cylinder = (w < width_inner) | (w > width_outer)

        projection = vg.Math.get_projection(self.coords, vnormal)
        mask_disk =  np.abs(projection) > height

        self.arr[mask_cylinder | mask_disk] = 0


# //////////////////////////////////////////////////////////////////////////////
