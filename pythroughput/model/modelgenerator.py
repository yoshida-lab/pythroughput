# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import random
import math
import pymatgen
from pymatgen.io.vasp.inputs import Poscar

"""
Model generator.
"""

logger = logging.getLogger(__name__)


class ModelGenerator(object):
    """
    Model generator. It can support symmetrical and unsymmetrical modification of models.
    In unsymmetrical modification, modifing cell length, cell shape (three angles),
    atom coordinates with varidating fractional coordinates and swapping species are containd.
    
    By treating models, we use the pymatgen.Structure class with dict-type representation.
    You can input any format of structure date it can be used in pymatgen package
    and treating it as pymatgen.Strucuture object in this class.
    Modified models can be output as pymatgen.Strucuture object.
    
    Arguments
    ---------
    struct: (fmt)
        The original structure.
    fmt: str
        The format of "struct".
    
    Parameters
    ----------
    struct: pymatgen.Structure
        The original structure.
    struct_dict: dict
        The (modified) structure with dict-type representation.
    """
    
    def __init__(self, struct, fmt="Structure"):
        """
        Arguments
        ---------
        struct: (fmt)
            The original structure.
        fmt: str
            The format of "struct".
        """
        if fmt is "Structure":
            self.struct = struct
        else:
            self.struct = pymatgen.Structure.from_str(open(struct).read(), fmt=fmt)
        self.struct_dict = self.struct.as_dict()
    
    def get_struct(self):
        """
        Gets modified structure.
        
        Returns
        -------
        pymatgen.Structure
            The modified structure.
        """
        return pymatgen.Structure.from_dict(self.struct_dict)
    
    def modify_symmetrical(self, min=-0.01, max=0.01):
        """
        Modifies given strucuture symmetrical,
        namely, only cell length are isometrically changed.
        
        Arguments
        ---------
        min: float
            Minimum change of cell length (%/100).
        max: float
            Maximum change of cell length (%/100).
        
        Parameters
        ----------
        var: float
            The change of cell length given by random.
        
        Returns
        -------
        self: ModelGenerator
        """
        var = random.uniform(1.00+min, 1.00+max)
        for i, row in enumerate(self.struct_dict["lattice"]["matrix"]):
            for j, column in enumerate(row):
                column *= var
                self.struct_dict["lattice"]["matrix"][i][j] = column
        return self
    
    def modify_unsymmetrical(self, modify_cell=True, modify_shape=False, modify_atom=False,
                             swap_atom=False, cell_min=-0.01, cell_max=0.01,
                             atom_modify_prob=10.0, atom_min=-0.01, atom_max=0.01,
                             atom_swap_num=10, atom_swap_restrict=True):
        """
        Modifies given structure unsymmetrical.
        
        Arguments
        ---------
        modify_cell: bool
            If modify cell length.
        modify_shape: bool
            If modify cell shape.
        modify_atom: bool
            If modify atom fractional coordinates.
        swap_atom: bool
            If swap two atoms.
        cell_min: float
            Minimum of changing cell (%/100).
        cell_max: float
            Maxmum of changing cell (%/100).
        atom_modify_prob: float
            Probability of changing atom coordinates (%).
        atom_min: float
            Minimum of changing atom fractional coordinates.
        atom_max: float
            Maxmum of changing atom fractional coordinates (%).
        atom_swap_num: int
            Number of swapping two atoms.
        atom_swap_resrtict: bool
            If swapping restrict to metal(non-metal)-metal(non-metal).
        
        Returns
        -------
        self
        """
        if modify_cell is True:
            self.modify_cell(modify_shape, cell_min, cell_max)
        if modify_atom is True:
            self.modify_atom(atom_modify_prob, atom_min, atom_max)
        if swap_atom is True:
            self.swap_atom(atom_swap_num, atom_swap_restrict)
        return self
    
    def modify_cell(self, modify_shape, min, max):
        """
        Modifies cell length and shape.
        
        Arguments
        ---------
        modify_shape: bool
            If modify cell shape.
        min: float
            Minimum of changing cell (%/100).
        max: float
            Maximum of changing cell (%/100).
        """
        for i, row in enumerate(self.struct_dict["lattice"]["matrix"]):
            if modify_shape is True:
                for j, column in enumerate(row):
                    # var_shape is to change cell shape if the elements of axis is 0.0.
                    # The validity of this type of changing shape is under testing.
                    var_shape = random.uniform(min, max)
                    var_cell = random.uniform(1.00+min, 1.00+max)
                    column += var_shape
                    column *= var_cell
                    self.struct_dict["lattice"]["matrix"][i][j] = column
            else:
                var = random.uniform(1.00+min, 1.00+max)
                for j, column in enumerate(row):
                    column *= var
                    self.struct_dict["lattice"]["matrix"][i][j] = column
    
    def modify_atom(self, prob, min, max):
        """
        Modify atomic coordinates.
        
        Arguments
        ---------
        prob: float
            Probability of changing atomic fractional coordinates (%).
        min: float
            Minimum of changing atomic fractional coordinates (%).
        max: float
            Maximum of changing atomic fractional coordinates (%).
        """
        for i, site in enumerate(self.struct_dict["sites"]):
            for j, coord in enumerate(site["abc"]):
                if random.uniform(0.0, 100.0) < prob:
                    coord += random.uniform(min, max)
                    self.struct_dict["sites"][i]["abc"][j] = coord
    
    def swap_atom(self, num, restrict):
        """
        Swaps to atoms.
        
        Arguments
        ---------
        num: int
            Number of swapping two atoms.
        restrict: bool
            If swapping restrict to metal(non-metal)-metal(non-metal).
        
        Parameters
        ----------
        atom1, atom2: int
            The index of two swapping atoms.
        dummy: pymatgen.Structure.element
            The intermidiate of swapping atoms.
        """
        for i in range(num):
            atom1, atom2 = self.select_atoms(restrict)
            dummy = self.struct_dict["sites"][atom1]["species"][0]["element"]
            self.struct_dict["sites"][atom1]["species"][0]["element"] = \
                self.struct_dict["sites"][atom2]["species"][0]["element"]
            self.struct_dict["sites"][atom2]["species"][0]["element"] = dummy
        self.sort_atom()
    
    def select_atoms(self, restrict):
        """
        Selects two swapping atoms.
        
        Arguments
        ---------
        restrict: bool
            If swapping restrict to metal(non-metal)-metal(non-metal).
        
        Parameters
        ----------
        atom1, atom2: int
            The index of two swapping atoms.
        
        Returns
        -------
        atom1, atom2: int
            The index of two swapping atoms.
        """
        while True:
            atom1 = random.randint(0, len(self.struct_dict["sites"])-1)
            atom2 = random.randint(0, len(self.struct_dict["sites"])-1)
            if not restrict or (self.is_metal(atom1) is self.is_metal(atom2)):
                break
        return atom1, atom2
    
    def is_metal(self, atom):
        """
        If given atom is metal.
        
        Arguments
        ---------
        atom: int
            The index of a swapping atoms.
        
        Returns
        -------
        bool
            If given atom is metal.
        """
        return self.struct_dict["sites"][atom]["species"][0]["element"] not in [
            "H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl",
            "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe", "At", "Rn"]
    
    def sort_atom(self):
        """
        Sorts atoms in struct_dict with key of atomic species.
        """
        self.struct_dict["sites"] = \
            sorted(self.struct_dict["sites"], key=lambda x:x["species"][0]["element"])
    
    def make_enough_large_supercell(self, num=10):
        """
        Makes enough large supercell for effective swapping.
        
        Arguments
        ---------
        num: int
            The lower limit of the number of atoms in the cell.
            If the number is lower than it, make shortest lattice vector
            twice larger by using pymatgen.Structure.make_supercell() method.
        
        Parameters
        ----------
        supercell: list
            The magnification rate of supercell,
            [2, 1, 1], [1, 2, 1] or [1, 1, 2].
        lattice: list
            The list of lattice length.
        min_lattice: int
            The index of lattice with minimum length.
        
        Returns
        -------
        bool
            If structure is largened.
        """
        if not self.struct.num_sites < num:
            return False
        while True:
            supercell = [1, 1, 1]
            lattice = [self.struct.lattice.a, self.struct.lattice.b, self.struct.lattice.c]
            min_lattice = lattice.index(min(lattice))
            supercell[min_lattice] = 2
            self.struct.make_supercell(supercell)
            if self.struct.num_sites < num:
                self.struct = self.get_struct()
                return True
    
    def reset_struct(self):
        """
        Resets modified structure to the initial one.
        
        Returns
        -------
        pymatgen.Structure
            Initial structure.
        """
        self.struct_dict = self.struct.as_dict()
        return self.struct
    
    def check_constrains(self, a_min=-1.0, a_max=100.0, b_min=-1.0, b_max=100.0,
                         c_min=-1.0, c_max=100.0, alpha_min=0.0, alpha_max=360.0,
                         beta_min=0.0, beta_max=360.0, gamma_min=0.0, gamma_max=360.0):
        """
        Checks constrains for modified structure.
        
        Arguments
        ---------
        a_min, a_max, b_min, b_max, c_min, c_max: float
            The constrains for cell length, a, b and c.
        alpha_min, alpha_max, beta_min, beta_max, gamma_min, gamma_max: float
            The constrains for cell angle alpha, beta and gamma.
        
        Returns
        -------
        bool
            If constrains holds.
        """
        return (a_min < self.struct_dict["lattice"]["a"] < a_max and
                b_min < self.struct_dict["lattice"]["b"] < b_max and
                c_min < self.struct_dict["lattice"]["c"] < c_max and
                alpha_min < self.struct_dict["lattice"]["alpha"] < alpha_max and
                beta_min  < self.struct_dict["lattice"]["beta"]  < beta_max  and
                gamma_min < self.struct_dict["lattice"]["gamma"] < gamma_max )
