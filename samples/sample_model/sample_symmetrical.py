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
The test for modify_symmetrical() method.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    gen = ModelGenerator("inputs/cif/Al2O3_hR30_R-3c_167.cif", fmt="cif")
    for i in range(10):
        while True:
            gen.modify_symmetrical(min=-0.05, max=0.05)
            constrains = gen.check_constrains(alpha_min=89.9, alpha_max=90.1,
                                              beta_min=89.9, beta_max=90.1,
                                              gamma_min=119.9, gamma_max=120.1)
            if constrains is True:
                break
        with open("outputs/symmetrical/POSCAR"+str(i).zfill(2), mode="w") as file:
            file.writelines(str(Poscar(gen.get_struct())))
        gen.reset_struct()
