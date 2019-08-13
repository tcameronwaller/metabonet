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

import metabonet.utility as utility

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
    path_candidacy = os.path.join(directory, "candidacy")
    path_metabolites = os.path.join(path_candidacy, "metabolites.pickle")
    path_reactions = os.path.join(path_candidacy, "reactions.pickle")
    # Read information from file.
    with open(path_metabolites, "rb") as file_source:
        metabolites_candidacy = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions_candidacy = pickle.load(file_source)
    # Specify directories and files.
    path_measurement = os.path.join(directory, "measurement")
    # Read information from file.
    study_one = read_source_study(
        study="study_one",
        path=path_measurement
    )
    study_two = read_source_study(
        study="study_two",
        path=path_measurement
    )
    study_three = read_source_study(
        study="study_three",
        path=path_measurement
    )
    study_four = read_source_study(
        study="study_four",
        path=path_measurement
    )
    study_five = read_source_study(
        study="study_five",
        path=path_measurement
    )
    # Compile and return information.
    return {
        "reactions_candidacy": reactions_candidacy,
        "metabolites_candidacy": metabolites_candidacy,
        "study_one": study_one,
        "study_two": study_two,
        "study_three": study_three,
        "study_four": study_four,
        "study_five": study_five
    }


def read_source_study(study=None, path=None):
    """
    Reads and organizes source information from file.

    arguments:
        study (str): name of study
        path (str): path to directory

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_study = os.path.join(path, (study + ".pickle"))
    # Read information from file.
    with open(path_study, "rb") as file_source:
        summary = pickle.load(file_source)
    # Return information.
    return summary


def match_measurements_to_candidate_metabolites(
    study_one=None,
    study_two=None,
    study_three=None,
    study_four=None,
    study_five=None,
    metabolites_candidacy=None
):
    """
    Matches measurements to candidate metabolites.

    arguments:
        study_one (list<dict<str>>): information about measurements for
            analytes
        study_two (list<dict<str>>): information about measurements for
            analytes
        study_three (list<dict<str>>): information about measurements for
            analytes
        study_four (list<dict<str>>): information about measurements for
            analytes
        study_five (list<dict<str>>): information about measurements for
            analytes
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites

    raises:

    returns:
        (dict<dict>): information about candidate metabolites

    """

    # Integrate information about measurements with information about candidate
    # metabolites.

    metabolites_measurements = {}
    for metabolite_candidacy in metabolites_candidacy.values():
        metabolite = metabolite_candidacy["metabolite"]
        measurements = determine_metabolite_measurements(
            metabolite=metabolite,
            study_one=study_one,
            study_two=study_two,
            study_three=study_three,
            study_four=study_four,
            study_five=study_five
        )
        metabolite_candidacy["measurements"] = measurements
        metabolites_measurements[metabolite_candidacy["identifier"]] = (
            metabolite_candidacy
        )
    return metabolites_measurements


def determine_metabolite_measurements(
    metabolite=None,
    study_one=None,
    study_two=None,
    study_three=None,
    study_four=None,
    study_five=None
):
    """
    Determines information about measurements for a metabolite.

    arguments:
        metabolite (str): identifier of a metabolite
        study_one (list<dict<str>>): information about measurements for
            analytes
        study_two (list<dict<str>>): information about measurements for
            analytes
        study_three (list<dict<str>>): information about measurements for
            analytes
        study_four (list<dict<str>>): information about measurements for
            analytes
        study_five (list<dict<str>>): information about measurements for
            analytes

    raises:

    returns:
        (dict<dict>): information about measurements for a candidate metabolite

    """

    # Determine measurements for the metabolite from each study.
    measurements_one = determine_metabolite_study_measurements(
        metabolite=metabolite,
        study=study_one
    )
    measurements_two = determine_metabolite_study_measurements(
        metabolite=metabolite,
        study=study_two
    )
    measurements_three = determine_metabolite_study_measurements(
        metabolite=metabolite,
        study=study_three
    )
    measurements_four = determine_metabolite_study_measurements(
        metabolite=metabolite,
        study=study_four
    )
    measurements_five = determine_metabolite_study_measurements(
        metabolite=metabolite,
        study=study_five
    )
    # Compile and return information.
    return {
        "study_one": measurements_one,
        "study_two": measurements_two,
        "study_three": measurements_three,
        "study_four": measurements_four,
        "study_five": measurements_five
    }


def determine_metabolite_study_measurements(
    metabolite=None,
    study=None
):
    """
    Determines information about measurements for a metabolite.

    arguments:
        metabolite (str): identifier of a metabolite
        study (list<dict<str>>): information about measurements for analytes

    raises:

    returns:
        (dict): information about measurements for a candidate metabolite

    """

    measurements = []
    for analyte in study:
        metabolites = analyte["references"]["metabolite"]
        if metabolite in metabolites:
            # Analyte matches metabolite.
            significance = analyte["significance"]
            if significance:
                significance_text = "True"
            else:
                significance_text = "False"
            record = {
                "fold": analyte["fold"],
                "fold_log": analyte["fold_log"],
                "p_value": analyte["p_value"],
                "p_value_log": analyte["p_value_log"],
                "significance": significance_text
            }
            measurements.append(record)
    # Determine whether any measurements match the metabolite.
    if len(measurements) < 1:
        # Include an empty measurement for the metabolite.
        record = {
            "fold": "",
            "fold_log": "",
            "p_value": "",
            "p_value_log": "",
            "significance": ""
        }
        measurements.append(record)
    return measurements[0]


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


def convert_metabolites_text(metabolites=None):
    """
    Converts information about metabolites to text format.

    arguments:
        metabolites (dict<dict>): information about metabolites

    returns:
        (list<dict>): information about metabolites

    raises:

    """

    records = []
    for metabolite in metabolites.values():
        measurements = metabolite["measurements"]
        study_one = measurements["study_one"]
        study_two = measurements["study_two"]
        study_three = measurements["study_three"]
        study_four = measurements["study_four"]
        study_five = measurements["study_five"]
        record = {
            "identifier": metabolite["identifier"],
            "name": metabolite["name"],
            "metabolite": metabolite["metabolite"],
            "compartment": metabolite["compartment"],
            "reactions_candidacy": ";".join(
                metabolite["reactions_candidacy"]
            ),
            "reactions_candidacy_count": (
                metabolite["reactions_candidacy_count"]
            ),
            "study_one_fold": study_one["fold"],
            "study_one_fold_log": study_one["fold_log"],
            "study_one_p_value": study_one["p_value"],
            "study_one_p_value_log": study_one["p_value_log"],
            "study_one_significance": study_one["significance"],
            "study_two_fold": study_two["fold"],
            "study_two_fold_log": study_two["fold_log"],
            "study_two_p_value": study_two["p_value"],
            "study_two_p_value_log": study_two["p_value_log"],
            "study_two_significance": study_two["significance"],
            "study_three_fold": study_three["fold"],
            "study_three_fold_log": study_three["fold_log"],
            "study_three_p_value": study_three["p_value"],
            "study_three_p_value_log": study_three["p_value_log"],
            "study_three_significance": study_three["significance"],
            "study_four_fold": study_four["fold"],
            "study_four_fold_log": study_four["fold_log"],
            "study_four_p_value": study_four["p_value"],
            "study_four_p_value_log": study_four["p_value_log"],
            "study_four_significance": study_four["significance"],
            "study_five_fold": study_five["fold"],
            "study_five_fold_log": study_five["fold_log"],
            "study_five_p_value": study_five["p_value"],
            "study_five_p_value_log": study_five["p_value_log"],
            "study_five_significance": study_five["significance"],
        }
        records.append(record)
    records.sort(
        key=lambda record: record["reactions_candidacy_count"],
        reverse=True
    )
    return records


def convert_reactions_text(reactions=None):
    """
    Converts information about reactions to text format.

    arguments:
        reactions (dict<dict>): information about reactions

    returns:
        (list<dict>): information about reactions

    raises:

    """

    records = []
    for reaction in reactions.values():
        record = {
            "identifier": reaction["identifier"],
            "reaction": reaction["reaction"],
            "name": reaction["name"],
            "replicates": ";".join(reaction["replicates"]),
            "metabolites_candidacy": ";".join(
                reaction["metabolites_candidacy"]
            ),
            "metabolites_candidacy_count": (
                reaction["metabolites_candidacy_count"]
            ),
            "measurement_one_fold": "",
            "measurement_one_log_fold": "",
            "measurement_one_p_value": "",
            "measurement_two_fold": "",
            "measurement_two_log_fold": "",
            "measurement_two_p_value": "",
            "measurement_three_fold": "",
            "measurement_three_log_fold": "",
            "measurement_three_p_value": "",
            "measurement_four_fold": "",
            "measurement_four_log_fold": "",
            "measurement_four_p_value": "",
            "measurement_five_fold": "",
            "measurement_five_log_fold": "",
            "measurement_five_p_value": ""
        }
        records.append(record)
    records.sort(
        key=lambda record: record["metabolites_candidacy_count"],
        reverse=True
    )
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
    path_network = os.path.join(directory, "network")
    path = os.path.join(path_network, "measurement")
    utility.confirm_path_directory(path)
    path_pickle = os.path.join(path, "metabolites.pickle")
    path_metabolites_text = os.path.join(path, "metabolites.tsv")
    path_reactions_text = os.path.join(path, "reactions.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["metabolites_measurement"], file_product)
    utility.write_file_table(
        information=information["metabolites_measurement_text"],
        path_file=path_metabolites_text,
        names=information["metabolites_measurement_text"][0].keys(),
        delimiter="\t"
    )
    if False:
        utility.write_file_table(
            information=information["reactions_measurement_text"],
            path_file=path_reactions_text,
            names=information["reactions_measurement_text"][0].keys(),
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
    metabolites_measurement = match_measurements_to_candidate_metabolites(
        study_one=source["study_one"],
        study_two=source["study_two"],
        study_three=source["study_three"],
        study_four=source["study_four"],
        study_five=source["study_five"],
        metabolites_candidacy=source["metabolites_candidacy"]
    )
    # Convert measurement information to table in text format.
    metabolites_measurement_text = convert_metabolites_text(
        metabolites=metabolites_measurement
    )
    reactions_measurement_text = convert_reactions_text(
        reactions=source["reactions_candidacy"]
    )
    # Compile information.
    information = {
        "metabolites_measurement": metabolites_measurement,
        "metabolites_measurement_text": metabolites_measurement_text,
        "reactions_measurement_text": reactions_measurement_text
    }
    #Write product information to file
    write_product(directory=directory, information=information)
