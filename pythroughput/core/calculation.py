# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import os
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
try:
    from pythroughput.core.calculation_vasp import Calculation_vasp
except ModuleNotFoundError:
    print("Raise ModuleNotFoundError: You cannot use VASP in this system!")
    pass

"""
Classes to performe high-throughput first-principles calculation.
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
                 input_path=None,
                 output_path=None,
                 **structs):
        """
        Arguments
        ---------
        calculator: dict
            Calculation configurations.
        input_path: str or None
            Path to input files other than structure files using in calculation.
            e.g.) path to potential files using in VASP calculation.
        output_path: str or None
            Path to output files. When it is None,
            output is omitted in GPAW calculation.
        structs: dict
            Dictionary of pymatgen.Structure object,
            which consists of the name of the structures
            as the keys and the atomic structures as the values.
        """
        self.structs = structs
        self.calculator = calculator
        self.input_path = input_path
        self.output_path = output_path
    
    def run(self, steps=1, package="gpaw"):
        """
        Runs high-throughput first-principles calculation.
        
        Arguments
        ---------
        steps: int
            Number of relaxation steps.
        package: str
            First-principles calculation package using in calculation.
        
        Returns
        -------
        results: dict
            Dictionary of calculation results, which consists of
            the name of structure as a key and the result of
            calculation as a value.
        """
        results = {}
        for struct_name, struct in self.structs.items():
            results[struct_name] = self._calc(struct_name, struct, steps, package)
        return self.results
    
    def _calc(self, struct_name, struct, steps, package):
        """
        Runs first-principles calculation.
        
        Arguments
        ---------
        struct_name: str
            Name of the structure.
        struct: pymatgen.Structure
            Atomic structure itself.
        steps:
            Number of relaxation steps.
        package: str
            First-principles calculation package using in calculation.
        
        Parameters
        ----------
        struct_calculator: dict
            Calculation configurations which are different
            depending on structures.
        
        Returns
        -------
        dict
            Calculation results return by get_results method.
        """
        struct_calculator = self._set_dafault_calculator(struct_name, struct)
        
        if package is "gpaw"
            try:
                return (
                    struct_name,
                    Calculation(struct_name, struct, struct_calculator).get_results(steps=steps)
                )
            except KohnShamConvergenceError:
                return (
                    struct_name,
                    {"results": "Unconverged"}
                )
        elif package is "vasp":
            return (
                struct_name,
                Calculation_vasp(struct_name, struct, struct_calculator, self.input_path).get_results(steps=steps)
        else:
            pass
    
    def _set_default_calculator(self, struct_name, struct):
        """
        Sets default calculation configulations which are different
        depending on structure, system size, and so on.
        
        Arguments
        ---------
        struct_name: str
            The name of the structure.
        struct: pymatgen.Structure
            Atomic structure itself.
        
        Returns
        -------
        struct_calculator: dict
            Calculation configurations which are different
            depending on structures.
        """
        struct_calculator = self.calculator
        if self.output_path is not None
            struct_calculator["txt"] = self.output_path + struct_name + ".txt"
        else:
            struct_calculator["txt"] = None
        
        # The optimal setting of kpts is under investigation.
        if struct.num_sites < 5:
            struct_calculator["kpts"] = {"size": (4, 4, 4)}
        elif struct.num_sites < 20:
            struct_calculator["kpts"] = {"size": (2, 2, 2)}
        else:
            struct_calculator["kpts"] = {"size": (1, 1, 1)}
        
        return struct_calculator
    

class Calculation(object):
    """
    Class to perform first-principles calculation with GPAW.
    """
    
    def __init__(self, struct_name, struct, calculator):
        """
        Arguments
        ---------
        struct_name:
            Name of the structure.
        struct: pymatgen.Structure
            Atomic structure itself.
        calculator: dict
            Calculation configurations.
        """
        self._struct_name = struct_name
        try:
            self._atom = AseAtomsAdaptor.get_atoms(struct, **struct.site_properties)
            self._atom.set_calculator(GPAW(**calculator))
        except ValueError:
            self._atom = None
    
    def get_results(self, steps=1):
        """
        Runs first-principles calculation with GPAW
        and gets calculation results.
        
        Arguments
        ---------
        steps: int
            Number of relaxation steps.
        
        Returns
        -------
        results: dict
            Dictionary of caluclation results.
        """
        if not self._is_not_invalid_struct():
            return self._struct_name, self._atom
        else:
            results = {}
            if not steps is 1:
                results["relax_struct"] = self._get_relax_struct(steps)
            else:
                pass
            results["total_energy"] = self._get_total_energy()
            return results
    
    def _is_not_invalid_struct(self):
        """
        Is not invalid structure.
        
        Returns
        -------
        bool
            True if self._atom is not None.
        """
        return self._atom is not None
    
    def _get_relax_struct(self, steps):
        """
        Gets relaxed structure by QuasiNewton calculation.
        
        Arguments
        ---------
        steps: int
            Number of relaxation steps.
        
        Returns
        -------
        pymatgen.Structure
            Relaxed structure.
        """
        QuasiNewton(self._atom, logfile=None).run(steps=steps)
        return AseAtomsAdaptor.get_structure(self._atom)
    
    def _get_total_energy(self):
        """
        Gets total energy (eV) of (relaxed) structure.
        
        Returns
        -------
        float
            Total energy (eV) of (relaxed) structure.
        """
        return self._atom.get_potential_energy()
