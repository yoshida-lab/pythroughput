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
from pythroughput.core.atomization import AtomizationCalculator

"""
Sample to perform high-throughput calculation of atomization energy.
"""

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    with open("outputs/csv/results.csv", mode="w") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(["struct_name", "total energy[eV]", "atomization energy[eV]"])
    
    structs = {}
    for struct_name in ("Fe4OF7",):
        structs[struct_name] = pymatgen.Structure.from_str(
            open("./inputs/" + struct_name + "/POSCAR").read(),
            fmt="poscar"
        )
    
    calc = PyHighThroughput(output_path="./outputs/calc/", **structs)
    calc.run()
    atomization = AtomizationCalculator(calc.results, ["O", "F", "Fe"])
    atomization.get_atomization_energy()
    
    with open("outputs/csv/results.csv", mode="a") as file:
        writer = csv.writer(file, lineterminator="\n")
        for struct_name, results in atomization.results.items():
            try:
                writer.writerow([
                    struct_name,
                    results["total_energy"],
                    results["atomization_energy"]
                ])
            except KeyError:
                writer.writerow([
                    struct_name,
                    "unconverged",
                    "undefined"
                ])
 
