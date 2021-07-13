import ipywidgets as ipw
from traitlets import Instance, Int, List, observe
from ase import Atoms, Atom
from mofchecker import MOFChecker
import tempfile
from itertools import chain


class MofCheckerWidget(ipw.VBox):
    """Widget that allows to check MOF structure correctness."""
    structure = Instance(Atoms, allow_none=True)
    selection = List(Int)

    def __init__(self, title=''):
        self.title = title
        button_check = ipw.Button(description="Analyse structure")
        button_check.on_click(self.check_structure)

        # Undercoodinated C.
        button_select_undercoord_c_ind = ipw.Button(description="Select these indices")
        button_select_undercoord_c_ind.on_click(self.select_underc_c_ind)
        self.l_undercoord_c_ind = ipw.HTML("<b>Undercoordinated C:</b> ")

        # Overcoordinated C.
        button_select_overcoord_c_ind = ipw.Button(description="Select these indices")
        button_select_overcoord_c_ind.on_click(self.select_overc_c_ind)
        self.l_overcoord_c_ind = ipw.HTML("<b>Overcoordinated C:</b> ")

        # Overcoordinated H.
        button_select_overcoord_h_ind = ipw.Button(description="Select these indices")
        button_select_overcoord_h_ind.on_click(self.select_overc_h_ind)
        self.l_overcoord_h_ind = ipw.HTML("<b>Overcoordinated H:</b> ")

        # Lone molecules.
        button_select_lone_mol_ind = ipw.Button(description="Select these indices")
        button_select_lone_mol_ind.on_click(self.select_lone_mol)
        self.l_lone_mol_ind = ipw.HTML("<b>Lone molecules:</b> ")

        # Add missing hydrogens.
        self.button_add_missing_hydrogens = ipw.Button(description="Add missing hydrogens", disabled=True)
        self.button_add_missing_hydrogens.on_click(self.add_missing_hydrogens)

        super().__init__(children=[button_check,
                                   ipw.HBox([button_select_overcoord_c_ind, self.l_overcoord_c_ind]),
                                   ipw.HBox([button_select_undercoord_c_ind, self.l_undercoord_c_ind]),
                                   ipw.HBox([button_select_overcoord_h_ind, self.l_overcoord_h_ind]),
                                   ipw.HBox([button_select_lone_mol_ind, self.l_lone_mol_ind]),
                                   self.button_add_missing_hydrogens,
                                  ])

    def check_structure(self, _=None):
        self.mfchk = MOFChecker.from_ase(self.structure, primitive=False)

        self.l_overcoord_c_ind.value = "<b>Overcoordinated C:</b> " + ', '.join(map(str, self.mfchk.overvalent_c_indices))
        self.l_undercoord_c_ind.value = "<b>Undercoordinated C:</b> " + ', '.join(map(str, self.mfchk.undercoordinated_c_indices))
        self.l_overcoord_h_ind.value = "<b>Overcoordinated H:</b> " + ', '.join(map(str, self.mfchk.overvalent_h_indices))
        self.l_lone_mol_ind.value = "<b>Lone molecules:</b> " + ', '.join(map(str, self.mfchk.lone_molecule_indices))
        self.button_add_missing_hydrogens.disabled = False

    def add_missing_hydrogens(self, _=None):
        """Add atoms."""
        atoms = self.structure.copy()
        selection = self.selection

        h_positions = self.mfchk.undercoordinated_c_candidate_positions + self.mfchk.undercoordinated_n_candidate_positions
        for pos in h_positions:
            atoms += Atom('H', pos)

        self.structure = atoms
        self.selection = selection

    def select_overc_c_ind(self, _=None):
        self.selection = self.mfchk.overvalent_c_indices

    def select_underc_c_ind(self, _=None):
        self.selection = self.mfchk.undercoordinated_c_indices

    def select_overc_h_ind(self, _=None):
        self.selection = self.mfchk.overvalent_h_indices

    def select_lone_mol(self, _=None):
        self.selection = list(chain.from_iterable(self.mfchk.lone_molecule_indices))

class CheckStructure(ipw.VBox):
    structure = Instance(Atoms, allow_none=True)

    def __init__(self, **kwargs):
        self.text = ipw.HTML()
        self.check_button = ipw.Button(description="Final Structure Check", button_style = 'warning')
        self.check_button.on_click(self.check_structure)
        super().__init__([self.check_button, self.text], **kwargs)

    def check_structure(self, _=None):
        if self.structure is None:
            self.text.value = ''
        else:
            self.text.value = ""
            self.mfchk = MOFChecker.from_ase(self.structure, primitive=False)
            issue_found = False
            for (check_name, check) in self.mfchk.checks.items():
                try:
                    if not check.is_ok:
                        issue_found = True
                    self.text.value += "\u2705 OK: " if check.is_ok else "\u274C Fail: "
                    self.text.value += f"{check.description} </br>"
                except:
                    pass
            self.check_button.button_style = 'danger' if issue_found else 'success'

    @observe('structure')
    def observe_structure(self, _=None):
        self.check_button.button_style = 'warning'
        self.text.value = ""

