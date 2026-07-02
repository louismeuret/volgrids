import warnings
import numpy as np

import volgrids as vg
import volgrids.smiffer as smf
import volgrids._vendors.molsimple as ms

from ._core.hbonds import SmifHBonds
from ._core.triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonors(SmifHBonds):
    MAX_BOND_DISTANCE = 1.5 # (Angstroms) maximum distance between donor and hydrogen to be considered bonded

    # --------------------------------------------------------------------------
    def __init__(self, mm: "smf.MoleculeManager"):
        super().__init__(mm)
        self.hbond_getter = smf.ParserChemTable.get_names_hbd
        self._kernel_hbd_free = vg.KernelGaussianBivariateAngleDist(
            radius = vg.CFG.param_hbd_free_dist_mu + vg.CFG.misc_kernel_gaussian_sigmas * vg.CFG.param_hbd_free_dist_sigma,
            deltas = self.mm.get_deltas(), dtype = vg.FLOAT_DTYPE, params = smf.PARAMS_HBD_FREE
        )
        self._kernel_hbd_fixed = vg.KernelGaussianBivariateAngleDist(
            radius = vg.CFG.param_hbd_fixed_dist_mu + vg.CFG.misc_kernel_gaussian_sigmas * vg.CFG.param_hbd_fixed_dist_sigma,
            deltas = self.mm.get_deltas(), dtype = vg.FLOAT_DTYPE, params = smf.PARAMS_HBD_FIXED
        )

        hydrogens = self.mm.atoms_all.select_non_water().select_hydrogens()
        if len(hydrogens) == 0:
            vg.CFG.smif_use_hydrogens = False

        self.bonds_to_h: dict[ms.Particle, ms.ParticleGroup] = self._attempt_to_guess_bonds(hydrogens)


    # --------------------------------------------------------------------------
    def can_be_interactor(self, triplet: Triplet) -> bool:
        if smf.ResnameStandard.is_prot(triplet.resname):
            if triplet.resname == "PRO": # donor only if there is no previous residue
                return not triplet.has_prev_res()

        if smf.ResnameStandard.is_nucleic(triplet.resname):
            if triplet.interactor == "O3'": # donor only if there is no next residue
                return not triplet.has_next_res()

            if triplet.interactor == "O5'": # donor only if there is no previous residue
                return not triplet.has_prev_res()

        return True # all other cases can be interactors


    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        if triplet.pos_head is not None: # head position is already set for succesful vg.CFG.smif_use_hydrogens iterations
            return

        triplet.set_pos_head()

        ############################### TAIL POSITION
        ### special cases for protein
        if smf.ResnameStandard.is_prot(triplet.resname):
            if triplet.interactor == "N": # tail points are in different residues
                if triplet.has_prev_res():
                    triplet.set_pos_tail_custom( # N of peptide bond
                        triplet.residue_prev,
                        triplet.residue_this,
                    )
                    self.kernel = self._kernel_hbd_fixed
                    return

                triplet.set_pos_tail( # N of N-terminus
                    tail_override = ("CA",) # usually tail is (C,CA) but here it needs to be overriden to only (CA,)
                )
                self.kernel = self._get_relevant_kernel(triplet)
                return

        triplet.set_pos_tail()


    # --------------------------------------------------------------------------
    def _iter_triplets(self):
        for triplet in super()._iter_triplets():
            if triplet.interactor in self.processed_interactors: continue

            if vg.CFG.smif_use_hydrogens:
                lst_atom_interactor = triplet.residue_this.select_name(triplet.interactor)
                if not lst_atom_interactor: continue
                atom_interactor = lst_atom_interactor[0] # the list should contain only one element, unless repeated atom names are present in the same residue

                for atom_h in self.bonds_to_h[atom_interactor]:
                    triplet.pos_tail = triplet.pos_interactor
                    triplet.pos_head = atom_h.get_position_numpy()
                    self.kernel = self._kernel_hbd_fixed
                    self.processed_interactors.add(triplet.interactor)
                    yield triplet

            if triplet.pos_head is None: # vg.CFG.smif_use_hydrogens falls back to "no-hydrogen" model if no hydrogens found
                self.kernel = self._get_relevant_kernel(triplet)
                yield triplet


    # --------------------------------------------------------------------------
    def _get_relevant_kernel(self, triplet: Triplet) -> vg.KernelGaussianBivariateAngleDist:
        return self._kernel_hbd_fixed if triplet.hbond_fixed else self._kernel_hbd_free


    # --------------------------------------------------------------------------
    def _attempt_to_guess_bonds(self, hydrogens: ms.ParticleGroup) -> list[ms.ParticleGroup]:
        if not vg.CFG.smif_use_hydrogens: return []

        coords_heavy = self.atoms.get_positions_numpy()
        coords_hydro = hydrogens.get_positions_numpy()
        mat_bonds = np.linalg.norm(
            coords_heavy[:, np.newaxis, :] - coords_hydro[np.newaxis, :, :],
            axis = 2
        ) < self.MAX_BOND_DISTANCE

        # return [hydrogens.select_mask(row) for row in mat_bonds]
        return {atom: hydrogens.select_mask(row) for atom,row in zip(self.atoms, mat_bonds)}


# //////////////////////////////////////////////////////////////////////////////
