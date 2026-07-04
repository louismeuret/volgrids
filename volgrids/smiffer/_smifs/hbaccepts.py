import volgrids as vg
import volgrids.smiffer as smf

from ._core.hbonds import SmifHBonds
from ._core.triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBAccepts(SmifHBonds):
    def __init__(self, mm: "smf.MoleculeManager"):
        super().__init__(mm)
        self.kernel = vg.KernelGaussianBivariateAngleDist(
            radius = vg.CFG.param_hba_dist_mu + vg.CFG.misc_kernel_gaussian_sigmas * vg.CFG.param_hba_dist_sigma,
            deltas = self.mm.get_deltas(), dtype = vg.FLOAT_DTYPE, params = smf.PARAMS_HBA
        )
        self.dict_triplets = mm.chemtable.names_hba


    # --------------------------------------------------------------------------
    def can_be_interactor(self, triplet: Triplet) -> bool:
        return True # acceptors can always be interactors, no special cases


    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        triplet.set_pos_head()

        ############################### TAIL POSITION
        ### special cases for RNA
        if smf.ResnameStandard.is_nucleic(triplet.resname):
            if triplet.interactor == "O3'": # tail points are in different residues
                triplet.set_pos_tail_custom(
                    triplet.residue_this,
                    triplet.residue_next,
                )
                return

        triplet.set_pos_tail()


# //////////////////////////////////////////////////////////////////////////////
