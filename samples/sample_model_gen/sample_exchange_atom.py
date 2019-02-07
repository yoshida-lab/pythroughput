# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the MIT License.

import logging
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import pythroughput
from pythroughput.model_gen.model_gen import Modelgen

"""
Sample for exchange_atom method of model_gen.py.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    modelgen_obj = Modelgen("./inputs/Al2O3.cif", format="cif")
    for i in range(3):
        modelgen_obj.exchange_atom(5)
        modelgen_obj.export_dict(filename="./outputs/sample_exchange_atom/POSCAR"+str(i+1))
