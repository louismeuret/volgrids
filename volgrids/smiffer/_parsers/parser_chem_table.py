from pathlib import Path

import volgrids as vg
import volgrids._vendors.molsimple as ms
import volgrids.smiffer as smf

# //////////////////////////////////////////////////////////////////////////////
class ParserChemTable:
    def __init__(self, path_table: Path):
        self.resnames: list[str] = []
        self.names_stk: dict[str, list[str]] = {}
        self.names_hba: dict[str, list[smf._smifs_core.Triplet]] = {}
        self.names_hbd: dict[str, list[smf._smifs_core.Triplet]] = {}

        self._atoms_hphob: dict[str, dict[str, float]] = {}
        self._parse_table(path_table)


    # --------------------------------------------------------------------------
    def get_atom_hphob(self, atom: ms.Particle) -> float | None:
        dict_resid = self._atoms_hphob.get(atom.resname)
        if dict_resid is None: return None
        return dict_resid.get(atom.name)


    # --------------------------------------------------------------------------
    def parse_atoms_hphobicity(self, data_ini: vg.ParserIni) -> None:
        for resname, str_groups in data_ini.iter_splitted_lines("HYDROPHOBICITY", sep = ':'):
            self._atoms_hphob[resname] = {}
            for group in str_groups.split():
                str_atoms, value = group.split('=')
                value = float(value)
                self._atoms_hphob[resname].update(**{
                    atom.strip() : value for atom in str_atoms.split(',')
                })


    # --------------------------------------------------------------------------
    def parse_names_stacking(self, data_ini: vg.ParserIni) -> None:
        for resname, str_cycles in data_ini.iter_splitted_lines("STACKING", sep = ':'):
            self.names_stk[resname] = [
                cycle.replace('-', ' ') for cycle in str_cycles.split()
            ]


    # --------------------------------------------------------------------------
    def parse_names_hbacceptors(self, data_ini: vg.ParserIni) -> None:
        for resname, str_triplets in data_ini.iter_splitted_lines("HBACCEPTORS", sep = ':'):
            triplets = list(map(self._parse_atoms_triplet, str_triplets.split()))
            self.names_hba[resname] = triplets


    # --------------------------------------------------------------------------
    def parse_names_hbdonors(self, data_ini: vg.ParserIni) -> None:
        for resname, str_triplets in data_ini.iter_splitted_lines("HBDONORS", sep = ':'):
            triplets = list(map(self._parse_atoms_triplet, str_triplets.split()))
            self.names_hbd[resname] = triplets


    # --------------------------------------------------------------------------
    def _parse_table(self, path_table: Path) -> None:
        """
        Populate the fields of the ParserChemTable instance by parsing the lines of the .chem table file.
        Should only be called once during initialization.
        """

        parser = vg.ParserIni.from_file(path_table)

        lst_resnames = parser.get("RESIDUE_NAMES")
        if not lst_resnames: raise ValueError("No selection query found in the table file.")
        self.resnames = lst_resnames[0].split()

        self.parse_atoms_hphobicity (parser)
        self.parse_names_stacking   (parser)
        self.parse_names_hbacceptors(parser)
        self.parse_names_hbdonors   (parser)


    # --------------------------------------------------------------------------
    @staticmethod
    def _parse_atoms_triplet(str_triplet: str) -> smf._smifs_core.Triplet:
        def _assert(condition: bool):
            assert condition, \
                f"Triplet '{str_triplet}' is not in the expected formats 'I=T->H' or 'I=T0.T1->H'."

        stripped = str_triplet.strip('!')
        hbond_fixed = stripped != str_triplet

        parts = stripped.split('=')
        _assert(len(parts) == 2)

        interactor, direction = parts
        parts = direction.split("->")
        _assert(len(parts) == 2)

        tail, head = parts
        _assert(tail and head and interactor)

        ### valid syntax?  | yes | no |
        ### 'i=t->h'       |  X  |    |
        ### 'i=t0.t1->h'   |  X  |    |
        ### 'i=t0.t1.t2->h'|  X  |    |
        ### 'i=t0.->h'     |     | X  |
        ### 'i=.t1->h'     |     | X  |
        tail_points = tail.split('.')
        _assert(len(tail_points) > 0 and all(tail_points))

        return smf._smifs_core.Triplet(interactor, tail_points, head, hbond_fixed)


# //////////////////////////////////////////////////////////////////////////////
