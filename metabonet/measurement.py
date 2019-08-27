"""
Author:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 4100, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of MetaboNet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports definition and analysis of custom metabolic networks.
    Copyright (C) 2019 Thomas Cameron Waller

    MetaboNet is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    MetaboNet is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with MetaboNet. If not, see <http://www.gnu.org/licenses/>.
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
    measurements = []
    path_measurement = os.path.join(directory, "measurement")
    if os.path.exists(path_measurement):
        for file in os.listdir(path_measurement):
            if file.endswith('.pickle'):
                path_study = os.path.join(path, file)
                # Read information from file.
                with open(path_study, "rb") as file_source:
                    summary = pickle.load(file_source)
                    measurements.append(summary)

    if len(measurements) < 1:
        measurements = None

    # Compile and return information.
    return {
        "reactions_candidacy": reactions_candidacy,
        "metabolites_candidacy": metabolites_candidacy,
        "studies": measurements
    }


def match_measurements_to_candidate_metabolites(
    studies=None,
    metabolites_candidacy=None
):
    """
    Matches measurements to candidate metabolites.

    arguments:
        studies (list<str>): information about measurements for
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
            studies=studies
        )
        metabolite_candidacy["measurements"] = measurements
        metabolites_measurements[metabolite_candidacy["identifier"]] = (
            metabolite_candidacy
        )
    return metabolites_measurements


def determine_metabolite_measurements(
    metabolite=None,
    studies=None
):
    """
    Determines information about measurements for a metabolite.

    arguments:
        metabolite (str): identifier of a metabolite
        studies (list<str>): information about measurements for
            analytes

    raises:

    returns:
        (dict<dict>): information about measurements for a candidate metabolite

    """

    # Determine measurements for the metabolite from each study.
    study_measurements = {}
    counter = 1

    for study_x in studies:

        measurement_x = determine_metabolite_study_measurements(
            metabolite=metabolite,
            study=study_x
        )
        study_name = str('study_') + str(counter)
        study_measurements[study_name] = measurement_x
        counter += 1

    # Compile and return information.
    return study_measurements


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

        }

        for key in measurements:
            study_x = measurements[key]

            key1 = str(key) + '_fold'
            record[key1] = study_x["fold"]

            key2 = str(key) + '_fold_log'
            record[key2] = study_x["fold_log"]

            key3 = str(key) + '_p_value'
            record[key3] = study_x["p_value"]

            key4 = str(key) + '_p_value_log'
            record[key4] = study_x["p_value_log"]

            key5 = str(key) + '_significance'
            record[key5] = study_x["significance"]

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
        studies=source["studies"],
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
