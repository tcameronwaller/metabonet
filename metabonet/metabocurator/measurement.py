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
    path_measurement = os.path.join(directory, "plasma_metabolites.tsv")
    path_extraction = os.path.join(directory, "extraction")
    path_hmdb = os.path.join(path_extraction, "hmdb_summary.pickle")
    path_conversion = os.path.join(directory, "conversion")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    # Read information from file.
    measurements = utility.read_file_table(
        path_file=path_measurement,
        names=[
            "set", "subset", "name", "reference_hmdb", "fold", "p_value",
            "q_value"
        ],
        delimiter="\t"
    )
    with open(path_hmdb, "rb") as file_source:
        summary_hmdb = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    # Compile and return information.
    return {
        "measurements": measurements,
        "summary_hmdb": summary_hmdb,
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
    for record in records[2:739]:
        name_original = record["name"]
        name_novel = name_original.replace("*", "")
        measurement = {
            "name": name_novel,
            "reference_hmdb": record["reference_hmdb"],
            "fold": float(record["fold"]),
            "p_value": float(record["p_value"])
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
    measurements_original=None, summary_hmdb=None
):
    """
    Enhances measurements' references to Human Metabolome Database (HMDB).

    arguments:
        measurements_original (list<dict>): information about measurements
        summary_hmdb (dict<dict>): information about metabolites from Human
            Metabolome Database (HMDB)

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        reference_hmdb = measurement["reference_hmdb"]
        name = measurement["name"]
        hmdb_keys = utility.match_hmdb_entries_by_identifiers_names(
            identifiers=[reference_hmdb],
            names=[name],
            summary_hmdb=summary_hmdb
        )
        del measurement["reference_hmdb"]
        measurement["hmdb"] = hmdb_keys
        measurements_novel.append(measurement)
    return measurements_novel


def match_measurements_to_metabolites(
    reference=None,
    measurements_original=None,
    metabolites=None
):
    """
    Matches measurements to metabolites.

    arguments:
        reference (str): name of attribute to use for match
        measurements_original (list<dict>): information about measurements
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        references_measurement = measurement[reference]
        # Find metabolites that match the reference.
        metabolites_matches = []
        for metabolite in metabolites.values():
            references_metabolite = metabolite["references"][reference]
            # Determine whether any measurement references match metabolite
            # references.
            matches = utility.filter_common_elements(
                list_one=references_measurement,
                list_two=references_metabolite
            )
            if len(matches) > 0:
                metabolites_matches.append(metabolite["identifier"])
        measurement["metabolites"] = metabolites_matches
        measurements_novel.append(measurement)
    return measurements_novel


def filter_measurements_metabolites(
    measurements_original=None
):
    """
    Filter measurements for those that map to metabolites.

    arguments:
        measurements_original (list<dict>): information about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        metabolites = measurement["metabolites"]
        if len(metabolites) > 0:
            measurements_novel.append(measurement)
    return measurements_novel


def filter_measurements_significance(
    p_value_threshold=None,
    measurements_original=None
):
    """
    Filter measurements by a threshold for significance.

    arguments:
        p_value_threshold (float): p value threshold for significance
        measurements_original (list<dict>): information about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    measurements_novel = []
    for measurement in measurements_original:
        p_value = measurement["p_value"]
        if p_value < p_value_threshold:
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
            "metabolites": ";".join(measurement["metabolites"])
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
    path_pickle = os.path.join(path, "measurements.pickle")
    path_text = os.path.join(path, "measurements.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["measurements"], file_product)
    utility.write_file_table(
        information=information["measurements_text"],
        path_file=path_text,
        names=information["measurements_text"][0].keys(),
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
    # Extract relevant information about measurements.
    measurements = extract_measurements(records=source["measurements"])
    # Match measurements to identifiers for Human Metabolome Database (HMDB).
    measurements_hmdb = enhance_measurements_hmdb_references(
        measurements_original=measurements,
        summary_hmdb=source["summary_hmdb"]
    )
    # Match measurements to metabolites.
    measurements_metabolites = match_measurements_to_metabolites(
        reference="hmdb",
        measurements_original=measurements_hmdb,
        metabolites=source["metabolites"]
    )
    # Filter measurements for those that map to metabolites.
    measurements_match = filter_measurements_metabolites(
        measurements_original=measurements_metabolites
    )
    # Calculate base-2 logarithm of fold change.
    measurements_log = calculate_measurements_log(
        measurements_original=measurements_match
    )
    # Filter analytes for those whose differences have p-values < 0.05.
    measurements_significance = filter_measurements_significance(
        p_value_threshold=0.05,
        measurements_original=measurements_log
    )
    # Convert measurement information to table in text format.
    measurements_text = convert_measurements_text(
        measurements=measurements_significance
    )
    # Compile information.
    information = {
        "measurements": measurements_significance,
        "measurements_text": measurements_text
    }
    #Write product information to file
    write_product(directory=directory, information=information)
