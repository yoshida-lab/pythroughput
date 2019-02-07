# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the MIT License.

import logging
import random
import math

import pymatgen
from pymatgen.io.vasp.inputs import Poscar

"""
Model generator in pythroughput module with pymatgen.
"""

logger = logging.getLogger(__name__)


class UndefinedFormatError(Exception):
    """
    Raises exception when user specified undefined format.
    """
    
    def __init__(self, format):
        self.format = format
    
    def __str__(self):
        return("Undefined format '{1}' is specified.".format(self.format))


class Modelgen(object):
    """
    Model generator class in pythroughput module with pymatgen.
    
    Parameters
    ----------
    struct: pymatgen.Structure
        Original structure data.
    struct_dict: dict
        Structure data which can be exchange atom positions.
    """
    
    def __init__(self, struct, format="structure",
                 filename="undefined", atom_num_limit=0,
                 output_original_struct_if_modified="True",
                 output_format="poscar", output_path="./"):
        """
        Arguments
        ---------
        struct: (Structure data)
            Original structure data, which can be pymatgen.Structure,
            poscar, cif, or any other what pymatgen can handle.
        format: str
            Format of read data.
        filename: str
            File name of read data.
        atom_num_limit: int
            The lower limit of the number of atoms in the unit cell.
        output_original_struct_if_modified: bool
            If this flag is True, output original structure
            if it is modified in make_supercell_enough_large method.
        output_format: str
            Format of output data of modified original structure.
        output_path: str
            Path to the output file.
        """
        if format is "structure":
            self.struct = struct
        else:
            self.struct = pymatgen.Structure.from_str(
                open(struct).read(),
                fmt=format
            )
        self.filename = filename
        
        struct_if_modified = self.make_supercell_enough_large(atom_num_limit)
        if output_original_struct_if_modified and struct_if_modified:
            if output_format is "poscar":
                with open(
                    output_path + "POSCAR_" + filename + "_modified",
                    mode="w"
                ) as file:
                    file.writelines(str(Poscar(self.struct)))
            else:
                raise UndefinedFormatError(format)
        
        self.struct_dict = self.struct.as_dict()
    
    def make_supercell_enough_large(self, atom_num_limit):
        """
        Makes supercell enough large to be our random exchanging effectively.
        
        Arguments
        ---------
        atom_num_limit: int
            The lower limit of the number of atoms in the unit cell.
            If the number is lower than it, make shortest lattice vector
            twice larger by using pymatgen.Structure.make_supercell method.
        """
        original_num_sites = self.struct.num_sites
        while self.struct.num_sites < atom_num_limit:
            lattice_len_dict = dict(zip(
                ["a", "b", "c"], 
                (self.struct.lattice.a, self.struct.lattice.b, self.struct.lattice.c)
            ))
            lattice_smallest_len_key = min(lattice_len_dict, key=lattice_len_dict.get)
            if lattice_smallest_len_key is "a":
                self.struct.make_supercell([2, 1, 1])
            elif lattice_smallest_len_key is "b":
                self.struct.make_supercell([1, 2, 1])
            else:
                self.struct.make_supercell([1, 1, 2])
        return not original_num_sites == self.struct.num_sites
    
    def exchange_atom(self, number=100, limit_exchange=False):
        """
        Run random exchange of atoms in structure.
        
        Arguments
        ---------
        number: int
            Number of exchange atoms (default: 100).
        restrict_exchange: bool
            If it is True, exchanging of atoms are limited in between
            metal-metal or non-metal-non-metal.
        """
        for i in range(number):
            while True:
                atom1 = random.randint(0, len(self.struct_dict["sites"])-1)
                atom2 = random.randint(0, len(self.struct_dict["sites"])-1)
                if not limit_exchange or self._is_metal(atom1) is self._is_metal(atom2):
                    break
            self._swap_atom(atom1, atom2)
        self._sort_atom()
        return self
    
    def _swap_atom(self, atom1, atom2):
        """
        Swap two atoms in struct_dict.
        
        Arguments
        ---------
        atom1, atom2: int
            Index of two swapping atoms in struct_dict.
        """
        element1 = self.struct_dict["sites"][atom1]["species"][0]["element"]
        element2 = self.struct_dict["sites"][atom2]["species"][0]["element"]
        self.struct_dict["sites"][atom1]["species"][0]["element"] = element2
        self.struct_dict["sites"][atom2]["species"][0]["element"] = element1
    
    def _sort_atom(self):
        """
        Sort atoms in struct_dict with key of atomic species.
        """
        sites_dict = self.struct_dict["sites"]
        sites_dict = sorted(
            sites_dict, key=lambda x:x["species"][0]["element"]
        )
        self.struct_dict["sites"] = sites_dict
    
    def _is_metal(self, atom):
        """
        Returns True if given atom is metal
        
        Arguments
        ---------
        atom: int
            Specie of atom in struct_dict.
        """
        return atom not in [
            "H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S",
            "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe", "At", "Rn"
        ]
    
    def modify_cell_size(self, probability=10.0, 
                         change_min=-1.0, change_max=1.0):
        """
        Modify cell size randomly.
        
        Arguments
        ---------
        probability: float
            Probability of changing cell size (%).
        change_min: float
            Minimum of changing cell size (% for the length of axes).
        change_max: float
            Maxmum of changing cell size (% for the length of axes).
        """
        len_axis = [] # 0: a-axis, 1: b-axis and 2: c-axis.
        for axis in self.struct_dict["lattice"]["matrix"]:
            len_axis.append(math.sqrt(axis[0]**2+axis[1]**2+axis[2]**2))
        
        for i, axis in enumerate(self.struct_dict["lattice"]["matrix"]):
            for j, element in enumerate(axis):
                if random.uniform(0.0, 100.0) < probability:
                    element += len_axis[i]*random.uniform(
                        change_min/100.0, change_max/100.0
                    )
                    self.struct_dict["lattice"]["matrix"][i][j] = element
                else:
                    pass
        return self
    
    def modify_atom_coord(self, probability=10.0,
                          change_min=-1.0, change_max=1.0):
        """
        Modify atomic coordinates randomly.
        
        Arguments
        ---------
        probability: float
            Probability of changing atomic coordinates (%).
        change_min: float
            Minimum of changing atomic coordinates (% for the length of axes).
        change_max: float
            Maxmum of changing atomic coordinates (% for the length of axes).
        """
        for i, site in enumerate(self.struct_dict["sites"]):
            for j, coord in enumerate(site["abc"]):
                if random.uniform(0.0, 100.0) < probability:
                    coord += random.uniform(
                        change_min/100.0, change_max/100.0
                    )
                    self.struct_dict["sites"][i]["abc"][j] = coord
                else:
                    pass
        return self
    
    def export_dict(self, format="poscar", filename="./POSCAR"):
        """
        Export edited structure.
        
        Arguments
        ---------
        fortmat: str
            Format of export data (default: "poscar").
        filename: str
            File name (directory) of exported file (default: "./POSCAR").
        """
        struct = pymatgen.Structure.from_dict(self.struct_dict)
        if format is "poscar":
            with open(filename, mode="w") as file:
                file.writelines(str(Poscar(struct)))
        else:
            raise UndefinedFormatError(format)
    

if __name__ == "__main__":
    pass
