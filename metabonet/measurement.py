"""
Extract information about metabolic sets and entities from MetaNetX.

Title:
    provision

Imports:
    os: This module from The Python Standard Library contains definitions of
        tools to interact with the operating system.
    sys: This module is from The Python Standard Library. It contains
        definitions of tools to interact with the interpreter.
    shutil: This module is from The Python Standard Library. It contains
        definitions of tools for file operations.
    importlib: This module is from The Python Standard library. It contains
        definitions of tools to import packages and modules.

Classes:
    This module does not contain any classes.

Exceptions:
    This module does not contain any exceptions.

Functions:
    ...

Author:
    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    Scientific Computing and Imaging Institute
    University Of Utah
    Room 4720 Warnock Engineering Building
    72 South Central Campus Drive
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of project Profondeur
    (https://github.com/tcameronwaller/profondeur/).

    Profondeur supports custom definition and visual exploration of metabolic
    networks.
    Copyright (C) 2018 Thomas Cameron Waller

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.
    If not, see <http://www.gnu.org/licenses/>.
"""


###############################################################################
# Notes

# The purpose of this procedure is to enhance information about metabolic sets
# and entities.

###############################################################################
# Installation and importation of packages and modules


# Packages and modules from the python standard library

import os
#import sys
import shutil
#import importlib
import csv
import copy
import pickle
import math

# Packages and modules from third parties

#import numpy
#import pandas
#import scipy

# Packages and modules from local source

import utility
import metabocurator.enhancement

###############################################################################
# Functionality


def read_source(directory=None):
    """
    Reads and organizes source information from file

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_measurement = os.path.join(directory, "measurement")
    path_measurements = os.path.join(path_measurement, "measurements.pickle")
    path_candidacy = os.path.join(directory, "candidacy")
    path_candidates = os.path.join(path_candidacy, "metabolites.pickle")
    # Read information from file.
    with open(path_measurements, "rb") as file_source:
        measurements = pickle.load(file_source)
    with open(path_candidates, "rb") as file_source:
        metabolites_candidacy = pickle.load(file_source)
    # Compile and return information.
    return {
        "measurements": measurements,
        "metabolites_candidacy": metabolites_candidacy
    }


def match_measurements_to_candidate_metabolites(
    measurements_original=None,
    metabolites_candidacy=None
):
    """
    Matches measurements to candidate metabolites.

    arguments:
        measurements_original (list<dict>): information about measurements
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        # Find candidate metabolites that match the measurement's metabolite.
        matches = []
        for candidate in metabolites_candidacy.values():
            # Determine whether candidate metabolite matches any of
            # measurement's metabolites.
            if candidate["metabolite"] in measurement["metabolites"]:
                matches.append(candidate["identifier"])
        measurement["candidates"] = matches
        measurements_novel.append(measurement)
    return measurements_novel


def filter_measurements_candidates(
    measurements_original=None
):
    """
    Filter measurements for those that map to candidate metabolites.

    arguments:
        measurements_original (list<dict>): information about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        candidates = measurement["candidates"]
        if len(candidates) > 0:
            measurements_novel.append(measurement)
    return measurements_novel


def convert_measurements_text(measurements=None):
    """
    Converts information about measurements to text format.

    arguments:
        measurements_original (list<dict>): information about measurements

    returns:
        (list<dict>): information about measurements

    raises:

    """

    records = []
    for measurement in measurements:
        record = {
            "name": measurement["name"],
            "fold": measurement["fold"],
            "log_fold": measurement["log_fold"],
            "p_value": measurement["p_value"],
            "hmdb": ";".join(measurement["hmdb"]),
            "metabolites": ";".join(measurement["metabolites"]),
            "candidates": ";".join(measurement["candidates"])
        }
        records.append(record)
    return records


def write_product(directory=None, information=None):
    """
    Writes product information to file

    arguments:
        directory (str): directory for product files
        information (object): information to write to file

    raises:

    returns:

    """

    # Specify directories and files.
    path = os.path.join(directory, "measurement")
    utility.confirm_path_directory(path)
    path_pickle = os.path.join(path, "measurements_candidacy.pickle")
    path_text = os.path.join(path, "measurements_candidacy_text.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["measurements_candidacy"], file_product)
    utility.write_file_table(
        information=information["measurements_candidacy_text"],
        path_file=path_text,
        names=information["measurements_candidacy_text"][0].keys(),
        delimiter="\t"
    )


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to extract relevant information from the
    Human Metabolome Database.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Match analytes to candidate metabolites.
    measurements_candidacy = match_measurements_to_candidate_metabolites(
        measurements_original=source["measurements"],
        metabolites_candidacy=source["metabolites_candidacy"]
    )
    # Filter measurements for those that map to candidate metabolites.
    measurements_match = filter_measurements_candidates(
        measurements_original=measurements_candidacy
    )
    # Convert measurement information to table in text format.
    measurements_candidacy_text = convert_measurements_text(
        measurements=measurements_match
    )
    # Compile information.
    information = {
        "measurements_candidacy": measurements_candidacy,
        "measurements_candidacy_text": measurements_candidacy_text
    }
    #Write product information to file
    write_product(directory=directory, information=information)
