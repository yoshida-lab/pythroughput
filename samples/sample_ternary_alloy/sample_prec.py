# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import csv
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import pymatgen
from pymatgen.io.vasp.inputs import Poscar
import pythroughput
from pythroughput.core.calculation import PyHighThroughput

"""
Sample to check precision and calculation time of GPAW calculation for ternary alloys.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    structs_path = {}
    
    for struct_name in ("Fe4OF7", "Fe4OF8", "Fe6OF11"):
        struct_path[struct_name] = ("./inputs/" + struct_name + "/POSCAR")
    
    structs = {}
    for struct_name, struct_path in structs_path.items():
        
        structs[struct_name] = pymatgen.Structure.from_str(
            open(struct_path + "_" + struct_name).read(),
            fmt="poscar"
        )
        
        structs[struct_name + "_var_1"] = pymatgen.Structure.from_str(
            open(struct_path + "1").read(),
            fmt="poscar"
        )
    
    pythroughput_obj = PyHighThroughput(
        txt_output=True,
        txt_output_path="./outputs/sample_prec/calc/",
        **structs
    )
    pythroughput_obj.run()
    
    with open("outputs/sample_prec/csv/results.csv", mode="w") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(["struct_name", "total energy [eV]"])
        
        for struct_name, results in pythroughput_obj.results.items():
            writer.writerow([struct_name, results["total_energy"]])
