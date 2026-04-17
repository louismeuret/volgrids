import io
import itertools
import logging
import warnings
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import numpy as np
import MDAnalysis as mda

import volgrids as vg
from volgrids.smiffer._parsers.parser_chem_table import ParserChemTable

# //////////////////////////////////////////////////////////////////////////////
class RNAResids:
    # --------------------------------------------------------------------------
    @classmethod
    def get_resids_bp_canonical(cls, path_pdb: Path) -> set[int]:
        structure2d = cls._get_rnapolis_struct2d(path_pdb)
        base_pairs = [bp for bp in structure2d.base_pairs if bp.saenger is not None]
        bp_full_names = [
            [bp.nt1.full_name, bp.nt2.full_name] for bp in base_pairs
            if bp.lw.value == "cWW" and bp.saenger.value in ("XIX", "XX")
        ]
        resids_0 = { int(bp[0][3:]) for bp in bp_full_names }
        resids_1 = { int(bp[1][3:]) for bp in bp_full_names }
        return resids_0 | resids_1


    # --------------------------------------------------------------------------
    @classmethod
    def get_resids_bp_wc(cls, path_pdb: str | Path) -> set[int]:
        """Detect Watson-Crick base pairs (AU, CG) using barnaba's H-bond counting logic.
        """
        _WC_MIN_HBONDS = {"AU": 2, "CG": 3}
        _HBOND_CUTOFF_SQ = 3.3 ** 2

        chem = ParserChemTable(vg.resolve_path_package("_tables/rna_simple_hb.chem"))
        u = cls._create_mda_universe_quiet(path_pdb)

        def _unique_names(getter, resname):
            return list(dict.fromkeys(t[0] for t in (getter(resname) or [])))

        def _positions(res_list, atom_names):
            # shape (n_res, n_atoms, 3); NaN where atom is absent
            arr = np.full((len(res_list), len(atom_names), 3), np.nan)
            for i, res in enumerate(res_list):
                for j, name in enumerate(atom_names):
                    sel = res.atoms.select_atoms(f"name {name}")
                    if len(sel):
                        arr[i, j] = sel[0].position
            return arr

        by_type = {}
        for base in ("A", "U", "C", "G"):
            res_list = [r for r in u.residues if r.resname.strip() == base]
            if not res_list:
                continue
            by_type[base] = {
                "res": res_list,
                "don": _positions(res_list, _unique_names(chem.get_names_hbd, base)),
                "acc": _positions(res_list, _unique_names(chem.get_names_hba, base)),
            }

        def _hbond_matrix(pos1, pos2):
            # pos1: (n1, k1, 3), pos2: (n2, k2, 3) → (n1, n2) count of close atom pairs
            # NaN distances propagate to NaN < cutoff == False, so missing atoms are ignored
            diff = pos1[:, np.newaxis, :, np.newaxis, :] - pos2[np.newaxis, :, np.newaxis, :, :]
            dist_sq = np.sum(diff ** 2, axis=-1)   # (n1, n2, k1, k2)
            return np.sum(dist_sq < _HBOND_CUTOFF_SQ, axis=(-2, -1))  # (n1, n2)

        wc_resids: set[int] = set()
        for b1, b2 in (("A", "U"), ("C", "G")):
            if b1 not in by_type or b2 not in by_type:
                continue
            n_hb = (
                _hbond_matrix(by_type[b1]["don"], by_type[b2]["acc"])
                + _hbond_matrix(by_type[b1]["acc"], by_type[b2]["don"])
            )  # (n_b1, n_b2)
            min_hb = _WC_MIN_HBONDS["".join(sorted([b1, b2]))]
            for i, j in np.argwhere(n_hb >= min_hb):
                wc_resids.add(int(by_type[b1]["res"][i].resid))
                wc_resids.add(int(by_type[b2]["res"][j].resid))

        return wc_resids


    # --------------------------------------------------------------------------
    @classmethod
    def get_all_resids(cls, path_pdb) -> set[int]:
        u = cls._create_mda_universe_quiet(path_pdb)
        return set(int(i) for i in u.residues.resids)


    # --------------------------------------------------------------------------
    @classmethod
    def get_resids_nonbp(cls, path_pdb: str | Path) -> str:
        path_pdb = Path(path_pdb)
        idxs_canonical = cls.get_resids_bp_canonical(path_pdb)
        idxs_all = cls.get_all_resids(path_pdb)
        idxs_nonpb = sorted(idxs_all - idxs_canonical)
        return ' '.join(str(i) for i in idxs_nonpb)


    # --------------------------------------------------------------------------
    @staticmethod
    def _create_mda_universe_quiet(path_pdb) -> mda.Universe:
        buf = io.StringIO()
        logger = logging.getLogger("MDAnalysis")
        old_level = logger.getEffectiveLevel()
        try:
            logger.setLevel(logging.ERROR)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with redirect_stdout(buf), redirect_stderr(buf):
                    u = mda.Universe(path_pdb)
        finally:
            logger.setLevel(old_level)
        return u


    # --------------------------------------------------------------------------
    @staticmethod
    def _get_rnapolis_struct2d(path_pdb: Path):
        from rnapolis.common import Structure2D
        from rnapolis.annotator import extract_base_interactions, handle_input_file, read_3d_structure
        file = handle_input_file(path_pdb)
        structure3d = read_3d_structure(file, None)
        base_interactions = extract_base_interactions(structure3d)
        structure2d, _ = structure3d.extract_secondary_structure(
            base_interactions, False, False
        )
        return structure2d


# //////////////////////////////////////////////////////////////////////////////
