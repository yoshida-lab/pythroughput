# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import sys
import os
import re
import pymatgen
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from pythroughput.core.calculation import PyHighThroughput

"""
Calculation module of atomization energy in pythroughput.
"""

logger = logging.getLogger(__name__)


class AtomizationCalculator(object):
    """
    Calculation module of atomization energy in pythroughput.
    
    This class gets csv which have total energy as an argument
    and calculate atomization energy by the total energy minus
    the energy of standard structures.
    
    The energy of standard structures are calculated in this class
    so you must give this class a calculator, which is dictionary
    in ase format, same as our module of pythroughput.
    """
    
    def __init__(self,
                 calculation_results,
                 species,
                 calculator={
                     "maxiter": 300,
                     "xc": "PBE",
                     "txt": None,
                     "kpts": {"size": (4, 4, 4)}
                 },
                 standard_struct_path=os.path.dirname(__file__)+"/standard_struct/"):
        """
        Arguments
        ---------
        calculation_result: dict
            Dictionary of calculation results, which are
            given by PyHighThroughput.results.
        species: list
            List of species.
        calculator: dict
            Calculation configurations.
        standard_struct_path: str
            Path to standard structures.
            Default: "./standard_struct/"
        """
        self.results = calculation_results
        self.species = species
        self.calculator = calculator
        self.standard_struct_path = standard_struct_path
    
    def get_atomization_energy(self, steps=1, package="gpaw", input_path=None):
        """
        Gets atomization energy.
        
        Arguments
        ---------
        package: str
            Calculation package using in calculation.
            Default: "gpaw"
        input_path: str or None
            Path to input files other than structure files using in calculation.
            e.g.) path to potential files using in VASP calculation.
        """
        self.standard_energy = self.calc_standard_energy(steps, package, input_path)
        self.calc_atomization_energy()
    
    def calc_standard_energy(self, steps, package, input_path=None):
        """
        Calculates standard energy.
        
        Arguments
        ---------
        package: str
            Calculation package using in calculation.
        
        Parameters
        ----------
        structs: list
            List of standard structures in pymatgen.Structure format.
        
        Returns
        -------
        list
            Total energies of standard structures you need.
        """
        if package is "vasp" and input_path is None:
            print("Error: In vasp calculation, you must set path to POTCARs as input_path.")
            return
        
        structs = {}
        
        for specie in self.species:
            structs[specie] = pymatgen.Structure.from_str(
                open(self.standard_struct_path+str(specie)+"/POSCAR").read(),
                fmt="poscar"
            )
        
        standard_energy = PyHighThroughput(
            input_path=input_path,
            output_path=None,
            **structs
        )
        standard_energy.run(steps=steps, package=package)
        return standard_energy.results
    
    def calc_atomization_energy(self):
        """
        Calculates atomization energy.
        """
        for struct_name, result in self.results.items():
            
            standard_energy_sum = 0
            
            try:
                for formula in result["formula"].split(" "):
                    specie = re.search("\D+", formula).group(0)
                    specie_num = re.search("\d+", formula).group(0)
                    standard_num = re.search(
                        "\d+", self.standard_energy[specie]["formula"]).group(0)
                    standard_energy_sum += (
                        self.standard_energy[specie]["total_energy"] / 
                        float(standard_num) * float(specie_num))
            
                if result.get("initial_energy") is not None:
                    result["initial_atomization"] = result["initial_energy"] - standard_energy_sum
                result["atomization_energy"] = result["total_energy"] - standard_energy_sum
            except KeyError:
                result["error"] = "KeyError"
