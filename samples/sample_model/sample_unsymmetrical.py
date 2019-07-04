# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import pymatgen
from pymatgen.io.vasp.inputs import Poscar
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model.modelgenerator import ModelGenerator

"""
The test for modify_unsymmetrical() method.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    gen = ModelGenerator("inputs/cif/Al2O3_hR30_R-3c_167.cif", fmt="cif")
    
    # cell length, only
    for i in range(3):
        gen.modify_unsymmetrical(cell_min=-0.05, cell_max=0.05)
        with open("outputs/unsymmetrical/POSCAR_len"+str(i).zfill(2), mode="w") as file:
            file.writelines(str(Poscar(gen.get_struct())))
        gen.reset_struct()
    
    # cell length and shape
    for i in range(3):
        gen.modify_unsymmetrical(modify_shape=True, cell_min=-0.05, cell_max=0.05)
        with open("outputs/unsymmetrical/POSCAR_shp"+str(i).zfill(2), mode="w") as file:
            file.writelines(str(Poscar(gen.get_struct())))
        gen.reset_struct()
    
    # atom coords, only
    for i in range(3):
        gen.modify_unsymmetrical(modify_cell=False, modify_atom=True,
                                 atom_min=-0.05, atom_max=0.05)
        with open("outputs/unsymmetrical/POSCAR_atm"+str(i).zfill(2), mode="w") as file:
            file.writelines(str(Poscar(gen.get_struct())))
        gen.reset_struct()
    
    # atom swap, only
    for i in range(3):
        gen.modify_unsymmetrical(modify_cell=False, swap_atom=True)
        with open("outputs/unsymmetrical/POSCAR_swp"+str(i).zfill(2), mode="w") as file:
            file.writelines(str(Poscar(gen.get_struct())))
        gen.reset_struct()
