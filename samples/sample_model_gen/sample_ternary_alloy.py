# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the MIT License.

import logging
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import pymatgen
import pythroughput
from pythroughput.model_gen.model_gen import Modelgen

"""
Practical sample for model_gen.py to generate models of ternary_alloys.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    structs = {}
    
    for struct_name in ("Fe4OF7", "Fe4OF8", "Fe6OF11"):
        structs[struct_name] = pymatgen.Structure.from_str(
            open("./inputs/POSCAR_" + struct_name).read(),
            fmt="poscar"
        )
    
    for struct_name, struct in structs.items():
        for i in range(33):
            modelgen_obj = Modelgen(
                struct,
                filename=struct_name,
                atom_num_limit=10,
                output_format="poscar",
                output_path="./outputs/sample_ternary_alloys/" + struct_name + "/"
            )
            
            modelgen_obj.modify_atom_coord()
            modelgen_obj.modify_cell_size()
            modelgen_obj.export_dict(
                filename="./outputs/sample_ternary_alloys/"+struct_name + "/POSCAR"+ str(i+1)
            )
