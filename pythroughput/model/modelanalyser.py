# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import numpy
import pymatgen
from pymatgen.io.vasp.inputs import Poscar

"""
Model analyser.
"""

logger = logging.getLogger(__name__)


class ModelAnalyser(object):
    """
    Model analyser.
    
    Arguments
    ---------
    structs: dict
        dict of analysed models, which format will be specified by "fmt".
    struct_stable: (fmt)
        Stable structure to analyse model, which format will be specified by "fmt".
    fmt: str
        structure format, default is pymatgen.Structure.
    
    Parameters
    ----------
    structs: dict
        dict of analysed models, format with pymatgen.Structure.
    struct_stable: pymatgen.Structure (or None)
        Stable structure to analyse model.
    """
    
    def __init__(self, structs, struct_stable=None, fmt="Structure"):
        if fmt == "Structure":
            self.structs = structs
            self.struct_stable = struct_stable
        else:
            self.structs = {}
            for struct_name, struct in structs.items():
                self.structs[struct_name] = pymatgen.Structure.from_str(
                    open(struct).read(), fmt=fmt
                )
            if struct_stable is not None:
                self.struct_stable = pymatgen.Structure.from_str(
                    open(struct_stable).read, fmt=fmt
                )
            else:
                self.struct_stable = struct_stable
    
    def calc_euclid_metric(self):
        """
        Analyse the metrics of models by calculationg euclidian distance.
        """
        metrics = {}
        for struct_name, struct in self.structs.items():
            metric = 0
            for i in range(3):
                metric += numpy.linalg.norm(
                    struct.lattice.matrix[i] - self.struct_stable.lattice.matrix[i]
                )
            for j in range(struct.num_sites):
                metric += numpy.linalg.norm(
                    numpy.array([struct.sites[j].a,
                                 struct.sites[j].b,
                                 struct.sites[j].c])
                    - numpy.array([self.struct_stable.sites[j].a,
                                   self.struct_stable.sites[j].b,
                                   self.struct_stable.sites[j].c])
                    )
            metrics[struct_name] = metric
        return metrics
