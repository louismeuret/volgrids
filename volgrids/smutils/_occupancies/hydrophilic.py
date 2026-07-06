import volgrids as vg
import volgrids.smiffer as smf

# //////////////////////////////////////////////////////////////////////////////
class OgHydrophilic(smf.SmifHydrophilic):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = vg.CFG.og_hphil_radius,
            deltas = self.mm.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid):
        """Populate grid with spherical accessibility regions."""
        grid.reset()
        for atom, logp_value in self.iter_particles():
            if logp_value > 0: continue
            self.kernel.stamp(grid, atom.get_position_numpy())


# //////////////////////////////////////////////////////////////////////////////
