# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

from .context import model
from model.modelgenerator import ModelGenerator
import random
import unittest
import logging
import pymatgen
from pymatgen.io.vasp.inputs import Poscar

"""
Test for modelgenerator.py
"""

logger = logging.getLogger(__name__)

# Generates model generator object.
gen = ModelGenerator("tests/inputs/cif/Al2O3_hR30_R-3c_167.cif", fmt="cif")

class ModelGeneratorTestSuite(unittest.TestCase):
    """
    Test for modelgenerator.py
    """
    
    def test_absolute_truth_and_meaning(self):
        """
        Absolute truth and meaning.
        """
        assert True
    
    def test_symmetric_model_generation(self):
        """
        Test for generating symmetric model.
        """
        random.seed(0)
        gen.modify_symmetrical(min=-0.01, max=0.01)
        struct = pymatgen.Structure.from_str(
            open("tests/samples/POSCAR_sym").read(), fmt="poscar"
        )
        assert str(Poscar(gen.get_struct())) == str(Poscar(struct))
    
    # Never passed, yet.
    def test_unsymmetric_model_generation(self):
        """
        Test for generating unsymmetric model.
        """
        random.seed(0)
        gen.modify_unsymmetrical(modify_shape=True, cell_min=-0.01, cell_max=0.01)
        struct = pymatgen.Structure.from_str(
            open("tests/samples/POSCAR_unsym").read(), fmt="poscar"
        )
        assert str(Poscar(gen.get_struct())) == str(Poscar(struct))
    

if __name__ == "__main__":
    unittest.main()
