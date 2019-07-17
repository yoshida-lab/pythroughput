# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

import logging
import csv

"""
CSV writer implimented to pythroughput package.
"""

logger = logging.getLogger(__name__)


class PyThroughCsv(object):
    """
    CSV writer implimented to pythroughput package.
    
    Arguments
    ---------
    filename: str
        CSV filename.
    results: list
        List of calculation results from pythroughput.Calculation classes.
    write_title: bool
        If title line is written automatically with formatting CSV file.
    parameters: list
            Spacification for writing results, if it is None,
            all the calculation results will be written.
    
    Parameters
    ----------
    filename: str
        CSV filename.
    results: list
        List of calculation results fron pythroughput.Calculation classes.
    """
    
    def __init__(self, filename, results, write_title=True, parameters=None):
        self.filename = filename
        self.results = results
        if parameters is None:
            self.parameters = list(list(self.results.values())[0].keys())
        else:
            self.parameters = parameters
        if write_title is True:
            with open(self.filename, mode="w") as file:
                writer = csv.writer(file, lineterminator="\n")
                writer.writerow(self.get_title_line())
    
    def get_title_line(self):
        """
        Gets first line of csvs.
        
        Returns
        -------
        title: line
            Title line.
        """
        title_line = []
        
        for parameter in self.parameters:
            title_line.append(parameter)
        
        return title_line
    
    def write_values(self):
        """
        Write values for csv file.
        """
        for result in self.results.values():
            with open(self.filename, mode="a") as file:
                writer = csv.writer(file, lineterminator="\n")
                writer.writerow(self.get_result_line(result))
    
    def get_result_line(self, result):
        """
        Gets result row to write for csvs.
        
        Arguments
        ---------
        result: dict
            Calculation results for each models.
        
        Returns
        -------
        result_line: list
            Result line.
        """
        result_line = []
        
        for parameter in self.parameters:
            try:
                result_line.append(result[parameter])
            except KeyError:
                result_line.append("Undefined key")
        
        return result_line
