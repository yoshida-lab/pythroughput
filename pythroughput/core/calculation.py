# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import os
# from joblib import Parallel, delayed
import pymatgen
from pymatgen.io.ase import AseAtomsAdaptor
try:
    from ase import Atom
    from ase.optimize import QuasiNewton
    from gpaw import GPAW
    from gpaw import KohnShamConvergenceError
except ModuleNotFoundError:
    print("Raise ModuleNotFoundError: You cannot use GPAW in this system!")
    pass

"""
Classes to performe high-throughput first-principles calculation with GPAW.
"""

logger = logging.getLogger(__name__)


class PyHighThroughput(object):
    """
    Class to performe high-throughput first-principles calculation with GPAW.
    
    Parameters
    ----------
    structs: pymatgen.Structure
        Calculating structures.
    results: dict
        dictionary of calculation results, which consist of
        physical properties as the keys and that value as the values.
    """
    
    def __init__(self,
                 calculator={
                     "maxiter": 300,
                     "xc": "PBE"
                 },
                 txt_output=False,
                 txt_output_path="./",
                 **structs):
        """
        Arguments
        ---------
        calculator: dict
            Calculation configuration.
        txt_output: bool
            Whether writing text output or not.
        txt_output_path: str
            Where to send text output if text_output is True.
        structs: dict
            dictionary of pymatgen.Structure object, which consists of
            the name of the structures as the keys and
            the atomic structures as the values.
        """
        self.structs = structs
        self.calculator = calculator
        self.txt_output = txt_output
        self.txt_output_path = txt_output_path
    
    def run(self, relaxation_steps=1):
        """
        Runs high-throughput first-principles calculation with GPAW.
        
        Arguments
        ---------
        relaxation_steps: int
            Number of steps of structural relaxation.
        """
        """Parallal computation is not supported in the current version.
        results = Parallel(n_jobs=1, verbose=0, timeout=None)([
            delayed(self._calc)(
                struct_name,
                struct,
                self.calculator,
                relaxation_steps,
                self.txt_output,
                self.txt_output_path
            ) for struct_name, struct in self.structs.items()
        ])
        """
        results = []
        for struct_name, struct in self.structs.items():
            results.append(self._calc(
                struct_name,
                struct,
                self.calculator,
                relaxation_steps,
                self.txt_output,
                self.txt_output_path
            ))
        self.results = dict(results)
        return self.results
    
    def _calc(self,
              struct_name,
              struct,
              calculator,
              steps,
              txt_output,
              txt_output_path):
        """
        Runs first-principles calculation with GPAW.
        
        Arguments
        ---------
        struct_name: str
            The name of the structure.
        struct: pymatgen.Structure
            Atomic structure.
        calculator: dict
            dictionary of calculation configuration.
        steps:
            Number of steps of structural relaxation.
        txt_output: bool
            Whether writing text output or not.
        txt_output_path: str
            Where to send text output if text_output is True.
        """
        calculator = self._update_calculator(
            struct_name,
            struct,
            calculator,
            txt_output,
            txt_output_path
        )
        try:
            return (
                struct_name,
                Calculation(struct_name, struct, calculator).get_results(steps)
            )
        except KohnShamConvergenceError:
            return (
                struct_name,
                {"results": "Unconverged"}
            )
        
    def _update_calculator(self,
                           struct_name,
                           struct,
                           calculator,
                           txt_output,
                           txt_output_path):
        """
        Updates calculator.
        
        Arguments
        ---------
        struct_name: str
            The name of the structure.
        struct: pymatgen.Structure
            Atomic structure.
        calculator: dict
            dictionary of calculation configuration.
        txt_output: bool
            Whether writing text output or not.
        txt_output_path: str
            Where to send text output if text_output is True.
        """
        if txt_output is True:
            calculator["txt"] = txt_output_path + struct_name + ".txt"
        else:
            calculator["txt"] = None

        if struct.num_sites < 5:
            calculator["kpts"] = {"size": (4, 4, 4)}
        elif struct.num_sites < 20:
            calculator["kpts"] = {"size": (2, 2, 2)}
        else:
            calculator["kpts"] = {"size": (1, 1, 1)}
        
        return calculator
    

class Calculation(object):
    """
    Class to perform first-principles calculation with GPAW.
    """
    
    def __init__(self, struct_name, struct, calculator):
        """
        Arguments
        ---------
        structs: pymatgen.Structure
            Structure data.
        calculator: dict
            Calculation config.
        """
        self._struct_name = struct_name
        try:
            self._atom = AseAtomsAdaptor.get_atoms(struct, **struct.site_properties)
            self._atom.set_calculator(GPAW(**calculator))
        except ValueError:
            self._atom = "Invalid structure"
    
    def get_results(self, steps=5):
        """
        Gets calculation results.
        
        Arguments
        ---------
        steps: int
            Number of steps of structural relaxation.
        
        Returns
        -------
        results: dict
            dictionary of caluclation results.
        """
        if not self._is_not_invalid_struct():
            return self._struct_name, self._atom
        else:
            results = {}
            if not steps is 1:
                results["relax_struct"] = self._get_relax_struct(steps)
            results["total_energy"] = self._get_total_energy()
            return results
    
    def _is_not_invalid_struct(self):
        """
        Return True if self._atom is not invalid structure
        """
        return self._atom is not "Invalid structure"
    
    def _get_relax_struct(self, steps):
        """
        Gets relaxed structure by QuasiNewton calculation.
        
        Returns
        -------
        pymatgen.Structure
            Relaxed structure.
        """
        QuasiNewton(self._atom, logfile=None).run(steps=steps)
        return AseAtomsAdaptor.get_structure(self._atom)
    
    def _get_total_energy(self):
        """
        Gets total energy (eV) of relaxed structure.
        
        Returns
        -------
        float
            Total energy (eV) of relaxed structure.
        """
        return self._atom.get_potential_energy()
