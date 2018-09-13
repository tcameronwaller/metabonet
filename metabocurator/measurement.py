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
import enhancement

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
    path_measurements = os.path.join(path_measurement, "plasma_metabolites.tsv")
    path_extrication = os.path.join(directory, "extrication")
    path_hmdb = os.path.join(path_extrication, "hmdb_summary.pickle")
    path_conversion = os.path.join(directory, "conversion")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    # Read information from file.
    measurements = utility.read_file_table(
        path_file=path_measurements,
        names=[
            "set", "subset", "name", "reference_hmdb", "fold", "p_value",
            "q_value"
        ],
        delimiter="\t"
    )
    with open(path_hmdb, "rb") as file_source:
        hmdb = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    # Compile and return information.
    return {
        "measurements": measurements,
        "hmdb": hmdb,
        "metabolites": metabolites
    }


def extract_measurements(records=None):
    """
    Extracts information about measurements.

    arguments:
        records (list<dict>): information from source about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements = []
    for record in records[2:]:
        measurement = {
            "name": record["name"],
            "reference_hmdb": record["reference_hmdb"],
            "fold": record["fold"],
            "p_value": record["p_value"]
        }
        measurements.append(measurement)
    return measurements


def calculate_measurements_log(measurements_original=None):
    """
    Calculates base-2 logarithm of fold change in measurements.

    arguments:
        records (list<dict>): information about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        fold = measurement["fold"]
        measurement["log_fold"] = math.log(fold, 2)
        measurements_novel.append(measurement)
    return measurements_novel


def enhance_measurements_hmdb_references(
    measurements_original=None, hmdb=None
):
    """
    Enhances measurements' references to Human Metabolome Database (HMDB).

    arguments:
        measurements_original (list<dict>): information about measurements
        hmdb (dict<dict>): information from HMDB

    raises:

    returns:
        (list<dict>): information about measurements

    """

    pass

def enhance_measurements_pubchem_references(
    measurements_original=None, hmdb=None
):
    """
    Enhances measurements' references to PubChem.

    arguments:
        measurements_original (list<dict>): information about measurements
        hmdb (dict<dict>): information from HMDB

    raises:

    returns:
        (list<dict>): information about measurements

    """

    pass


def enhance_measurements_metanetx_references(
    measurements_original=None, metabolites=None
):
    """
    Enhances measurements' references to MetaNetX.

    arguments:
        measurements_original (list<dict>): information about measurements
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (list<dict>): information about measurements

    """

    pass


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
    path = os.path.join(directory, "provision")
    utility.confirm_path_directory(path)
    path_pickle = os.path.join(path, "hmdb_summary.pickle")
    path_text = os.path.join(path, "hmdb_summary.tsv")
    path_midas = os.path.join(path, "midas_novel.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["summary_object"], file_product)
    utility.write_file_table(
        information=information["summary_list"],
        path_file=path_text,
        names=information["summary_list"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["midas_novel"],
        path_file=path_midas,
        names=information["midas_novel"][0].keys(),
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
    # TODO: read in the raw measurement table from the publication
    # TODO: read in the summary of hmdb
    # TODO: read in curated information about metabolites
    source = read_source(directory=directory)
    # Extract relevant information about measurements.
    measurements = extract_measurements(records=source["measurements"])
    # Calculate base-2 logarithm of fold change.
    measurements_log = calculate_measurements_log(
        measurements_original=measurements
    )
    # Match analytes to identifiers for Human Metabolome Database (HMDB),
    # PubChem, and MetaNetX.
    # TODO: match records by alternative HMDB identifiers (include these in HMDB summary)
    # TODO: match records by comparing name to synonymous names for HMDB entry (include these in HMDB summary)
    measurements_hmdb = enhance_measurements_hmdb_references(
        measurements_original=measurements_log,
        hmdb=source["hmdb"]
    )
    measurements_pubchem = enhance_measurements_pubchem_references(
        measurements_original=measurements_hmdb,
        hmdb=source["hmdb"]
    )
    measurements_metanetx = enhance_measurements_metanetx_references(
        measurements_original=measurements_pubchem,
        metabolites=source["metabolites"]
    )
    # Filter analytes for those whose differences have p-values < 0.05.

    # Convert measurement information to table in text format.


    # Compile information.
    information = {
        "midas_novel": midas_novel
    }
    #Write product information to file
    write_product(directory=directory, information=information)
