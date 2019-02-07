# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import csv
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import pymatgen
from pymatgen.io.vasp.inputs import Poscar
import pythroughput
from pythroughput.core.calculation import PyHighThroughput

"""
Sample to perform high-throughput calculation for basic metal structrues.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    structs = {}
    
    # read input structures for dictionary "structs".
    for filename in ("Al_fcc", "Mg_bcc", "Mg_hcp",
                     "Ti_bcc", "Ti_hcp", "NaCl", "CsCl"):
        
        structs[filename] = (
            pymatgen.Structure.from_str(
                open("inputs/" + filename + ".cif").read(),
                fmt="cif"
            )
        )
    
    pythroughput_obj = PyHighThroughput(**structs)
    pythroughput_obj.run()
    
    # output calculation results.
    with open("outputs/results.csv", mode="w") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(["struct_name", "total energy [eV]"])
        
        for struct_name, results in pythroughput_obj.results.items():
            writer.writerow([struct_name, results["total_energy"]])
