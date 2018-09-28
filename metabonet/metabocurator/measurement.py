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
import textwrap
import statistics

# Packages and modules from third parties

#import numpy
#import pandas

import scipy.stats

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
    # Read information from file.
    # Compile information.
    reference = read_source_reference(directory=directory)
    study_one = read_source_study_one(directory=directory)
    study_two = read_source_study_two(directory=directory)
    study_three_four = read_source_study_three_four(directory=directory)
    # Compile and return information.
    return {
        "reference": reference,
        "study_one": study_one,
        "study_two": study_two,
        "study_three_four": study_three_four
    }


def read_source_reference(directory=None):
    """
    Reads and organizes source information from file.

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_extraction = os.path.join(directory, "extraction")
    path_hmdb = os.path.join(path_extraction, "hmdb_summary.pickle")
    path_conversion = os.path.join(directory, "conversion")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    # Read information from file.
    with open(path_hmdb, "rb") as file_source:
        hmdb = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    # Compile and return information.
    return {
        "hmdb": hmdb,
        "metabolites": metabolites
    }


def read_source_study_zero(directory=None):
    """
    Reads and organizes source information from file.

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_measurement = os.path.join(directory, "metabolomic_measurements")
    path_study_zero = os.path.join(
        path_measurement, "karl_physiological-reports_2017"
    )
    path_measurements = os.path.join(
        path_study_zero, "measurements.tsv"
    )
    # Read information from file.
    measurements = utility.read_file_table(
        path_file=path_measurements,
        names=["set", "subset", "name", "hmdb", "fold", "p_value", "q_value"],
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "measurements": measurements
    }


def read_source_study_one(directory=None):
    """
    Reads and organizes source information from file.

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_measurement = os.path.join(directory, "metabolomic_measurements")
    path_study = os.path.join(
        path_measurement, "metabolomics-workbench_pr000058_st000061"
    )
    path_samples = os.path.join(path_study, "samples.tsv")
    path_analytes = os.path.join(path_study, "analytes.tsv")
    path_measurements = os.path.join(path_study, "measurements.tsv")
    # Read information from file.
    samples = utility.read_file_table(
        path_file=path_samples,
        names=None,
        delimiter="\t"
    )
    analytes = utility.read_file_table(
        path_file=path_analytes,
        names=None,
        delimiter="\t"
    )
    measurements = utility.read_file_table(
        path_file=path_measurements,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "samples": samples,
        "analytes": analytes,
        "measurements": measurements
    }


def read_source_study_two(directory=None):
    """
    Reads and organizes source information from file.

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_measurement = os.path.join(directory, "metabolomic_measurements")
    path_study = os.path.join(
        path_measurement, "metabolomics-workbench_pr000305_st000390"
    )
    path_samples = os.path.join(path_study, "samples.tsv")
    path_analytes = os.path.join(path_study, "analytes.tsv")
    path_measurements = os.path.join(path_study, "measurements.tsv")
    # Read information from file.
    samples = utility.read_file_table(
        path_file=path_samples,
        names=None,
        delimiter="\t"
    )
    analytes = utility.read_file_table(
        path_file=path_analytes,
        names=None,
        delimiter="\t"
    )
    measurements = utility.read_file_table(
        path_file=path_measurements,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "samples": samples,
        "analytes": analytes,
        "measurements": measurements
    }


def read_source_study_three_four(directory=None):
    """
    Reads and organizes source information from file.

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_measurement = os.path.join(directory, "metabolomic_measurements")
    path_study = os.path.join(
        path_measurement, "metabolomics-workbench_pr000322_st000412"
    )
    path_samples = os.path.join(path_study, "samples.tsv")
    path_analytes = os.path.join(path_study, "analytes.tsv")
    path_measurements = os.path.join(path_study, "measurements.tsv")
    # Read information from file.
    samples = utility.read_file_table(
        path_file=path_samples,
        names=None,
        delimiter="\t"
    )
    analytes = utility.read_file_table(
        path_file=path_analytes,
        names=None,
        delimiter="\t"
    )
    measurements = utility.read_file_table(
        path_file=path_measurements,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "samples": samples,
        "analytes": analytes,
        "measurements": measurements
    }


# Study zero.


def curate_measurements_study_zero(measurements=None):
    """
    Extracts information about measurements.

    arguments:
        measurements (list<dict>): information from source about measurements

    raises:

    returns:
        (list<dict>): information about measurements

    """

    # Extract relevant information about measurements.
    measurements = extract_measurements_study_one(
        records=source["measurements"]
    )
    # Match measurements to identifiers for Human Metabolome Database (HMDB).
    measurements_hmdb = enhance_measurements_hmdb_references(
        measurements_original=copy.deepcopy(measurements),
        summary_hmdb=source["summary_hmdb"]
    )
    # Match measurements to metabolites.
    measurements_metabolites = match_measurements_to_metabolites(
        reference="hmdb",
        measurements_original=copy.deepcopy(measurements_hmdb),
        metabolites=source["metabolites"]
    )
    # Filter measurements for those that map to metabolites.
    measurements_match = filter_measurements_metabolites(
        measurements_original=copy.deepcopy(measurements_metabolites)
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
    # Compile and return information.
    return {
        "measurements": measurements_significance,
        "measurements_text": measurements_text
    }


def extract_measurements_study_zero(records=None):
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
            "hmdb": record["hmdb"],
            "fold": float(record["fold"]),
            "p_value": float(record["p_value"])
        }
        measurements.append(measurement)
    return measurements


# Study with dependent pairs of samples.


def curate_measurements_study_pairs(
    group_numerator=None,
    group_denominator=None,
    samples=None,
    analytes=None,
    measurements=None,
    hmdb=None,
    metabolites=None
):
    """
    Curates information about metabolomic measurements from a study.

    arguments:
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        samples (list<dict<str>>): information about samples from a study
        analytes (list<dict<str>>): information about analytes from a study
        measurements (list<dict<str>>): information about measurements from a
            study
        hmdb (dict<dict>): information about metabolites from Human Metabolome
            Database (HMDB)
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (list<dict>): information about measurements

    """

    # Derive summary from analytes.
    summary = extract_analytes_summary(
        analytes=analytes
    )
    # Enhance analyte references.
    summary_reference = enhance_analytes_references(
        summary=summary,
        hmdb=hmdb
    )
    # Match analytes to metabolites.
    # Match by PubChem identifiers.
    summary_metabolite = match_analytes_to_metabolites(
        reference="pubchem",
        summary=summary_reference,
        metabolites=metabolites
    )
    # Determine fold changes.
    summary_fold = calculate_analytes_pairs_folds(
        summary=summary_metabolite,
        group_numerator=group_numerator,
        group_denominator=group_denominator,
        samples=samples,
        measurements=measurements
    )
    # Determine logarithms-base-2 of fold changes.
    summary_log = calculate_folds_logarithms(
        records=summary_fold
    )
    # Determine p-values.
    # Compare pairs of samples in both groups.
    # Apply pair t-test for dependent sample populations.
    summary_p = calculate_analytes_pairs_p_values(
        summary=summary_log,
        group_one=group_numerator,
        group_two=group_denominator,
        samples=samples,
        measurements=measurements
    )
    # Filter for anlaytes that match metabolites.
    summary_match = filter_analytes_metabolites(
        summary=copy.deepcopy(summary_p)
    )
    # Convert measurement information to table in text format.
    summary_text = convert_summary_text(summary=summary_match)
    # Report.
    print("analytes, measurements after curation...")
    report = prepare_curation_report(summary=summary_p)
    print(report)
    # Compile and return information.
    return {
        "summary": summary_match,
        "summary_text": summary_text
    }


# Study with independent samples.


def curate_measurements_study(
    group_numerator=None,
    group_denominator=None,
    samples=None,
    analytes=None,
    measurements=None,
    hmdb=None,
    metabolites=None
):
    """
    Curates information about metabolomic measurements from a study.

    arguments:
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        samples (list<dict<str>>): information about samples from a study
        analytes (list<dict<str>>): information about analytes from a study
        measurements (list<dict<str>>): information about measurements from a
            study
        hmdb (dict<dict>): information about metabolites from Human Metabolome
            Database (HMDB)
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (list<dict>): information about measurements

    """

    # Derive summary from analytes.
    summary = extract_analytes_summary(
        analytes=analytes
    )
    # Enhance analyte references.
    summary_reference = enhance_analytes_references(
        summary=summary,
        hmdb=hmdb
    )
    # Match analytes to metabolites.
    # Match by PubChem identifiers.
    summary_metabolite = match_analytes_to_metabolites(
        reference="pubchem",
        summary=summary_reference,
        metabolites=metabolites
    )
    # Determine fold changes.
    summary_fold = calculate_analytes_folds(
        summary=summary_metabolite,
        group_numerator=group_numerator,
        group_denominator=group_denominator,
        samples=samples,
        measurements=measurements
    )
    # Determine logarithms-base-2 of fold changes.
    summary_log = calculate_folds_logarithms(
        records=summary_fold
    )
    # Determine p-values.
    summary_p = calculate_analytes_p_values(
        summary=summary_log,
        group_one=group_numerator,
        group_two=group_denominator,
        samples=samples,
        measurements=measurements
    )
    # Filter for anlaytes that match metabolites.
    summary_match = filter_analytes_metabolites(
        summary=copy.deepcopy(summary_p)
    )
    # Convert measurement information to table in text format.
    summary_text = convert_summary_text(summary=summary_match)
    # Report.
    print("analytes, measurements after curation...")
    report = prepare_curation_report(summary=summary_p)
    print(report)
    # Compile and return information.
    return {
        "summary": summary_match,
        "summary_text": summary_text
    }


# General utility.


def extract_analytes_summary(analytes=None):
    """
    Extracts information about analytes.

    arguments:
        analytes (list<dict<str>>): information about analytes from a study

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    summary = []
    for analyte in analytes:
        identifier = analyte["identifier"]
        name = analyte["name"]
        metabolomics_workbench = analyte["reference_metabolomics_workbench"]
        if analyte["reference_pubchem"] == "-":
            pubchem = ""
        else:
            pubchem = analyte["reference_pubchem"]
        references = {
            "metabolomics_workbench": metabolomics_workbench,
            "pubchem": pubchem
        }
        record = {
            "identifier": identifier,
            "name": name,
            "references": references
        }
        summary.append(record)
    return summary


def enhance_analytes_references(
    summary=None, hmdb=None
):
    """
    Enhances analytes' references to Human Metabolome Database (HMDB) and to
    PubChem.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes
        hmdb (dict<dict>): information about metabolites from Human Metabolome
            Database (HMDB)

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    summary_novel = []
    for record in summary:
        name_one = record["identifier"]
        name_two = record["name"]
        hmdb_keys = utility.match_hmdb_entries_by_identifiers_names(
            identifiers=[],
            names=[name_one, name_two],
            summary_hmdb=hmdb
        )
        # Include references to HMDB.
        record["references"]["hmdb"] = hmdb_keys
        # Enhance references to PubChem.
        pubchem = []
        for key in hmdb_keys:
            hmdb_entry = hmdb[key]
            hmdb_pubchem = hmdb_entry["reference_pubchem"]
            pubchem.append(hmdb_pubchem)
        if len(record["references"]["pubchem"]) > 0:
            pubchem.append(record["references"]["pubchem"])
        pubchem_unique = utility.collect_unique_elements(pubchem)
        record["references"]["pubchem"] = pubchem_unique
        summary_novel.append(record)
    return summary_novel


def match_analytes_to_metabolites(
    reference=None,
    summary=None,
    metabolites=None
):
    """
    Matches measurements to metabolites.

    arguments:
        reference (str): name of attribute to use for match
        summary (list<dict<str>>): information about measurements for analytes
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    summary_novel = []
    for record in summary:
        references_record = [record["references"][reference]]
        # Find metabolites that match the record's reference.
        metabolites_matches = []
        for metabolite in metabolites.values():
            references_metabolite = metabolite["references"][reference]
            # Determine whether any references match.
            matches = utility.filter_common_elements(
                list_one=references_record,
                list_two=references_metabolite
            )
            if len(matches) > 0:
                metabolites_matches.append(metabolite["identifier"])
        record["references"]["metabolite"] = metabolites_matches
        summary_novel.append(record)
    return summary_novel


def calculate_analytes_folds(
    summary=None,
    group_numerator=None,
    group_denominator=None,
    samples=None,
    measurements=None
):
    """
    Determines the mean fold changes for each analyte between experimental
    groups.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        samples (list<dict<str>>): information about samples from a study
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    # Determine samples in each group.
    groups_samples = determine_groups_samples(
        group_one=group_numerator,
        group_two=group_denominator,
        samples=samples
    )
    # Determine mean fold changes for each analyte.
    summary_novel = []
    for record in summary:
        identifier = record["identifier"]
        fold = calculate_analyte_fold(
            identifier=identifier,
            group_numerator=group_numerator,
            group_denominator=group_denominator,
            groups_samples=groups_samples,
            measurements=measurements
        )
        record["fold"] = fold
        summary_novel.append(record)
    return summary_novel


def determine_groups_samples(
    group_one=None,
    group_two=None,
    samples=None
):
    """
    Determines the samples for each group for the same patient.

    arguments:
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        samples (list<dict<str>>): information about samples from a study

    raises:

    returns:
        (dict): samples from both groups

    """

    groups = [group_one, group_two]
    groups_samples = {
        group_one: [],
        group_two: []
    }
    for sample in samples:
        identifier = sample["identifier"]
        group = sample["group"]
        # Determine whether the sample belongs to a relevant group.
        if group in groups:
            groups_samples[group].append(identifier)
    return groups_samples


def calculate_analyte_fold(
    identifier=None,
    group_numerator=None,
    group_denominator=None,
    groups_samples=None,
    measurements=None
):
    """
    Determines the mean fold change between experimental groups for an analyte.

    arguments:
        identifier (str): identifier of an analyte
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        groups_samples (dict<list<str>>): samples from both groups
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (float): fold change

    """

    # Find measurements for analyte.
    measurements_analyte = find_analyte_measurements(
        identifier=identifier,
        measurements=measurements
    )
    # Determine fold changes for analyte's measurements.
    numerator_values = []
    for sample in groups_samples[group_numerator]:
        if (sample in measurements_analyte.keys()):
            value = float(measurements_analyte[sample])
            numerator_values.append(value)
    denominator_values = []
    for sample in groups_samples[group_denominator]:
        if (sample in measurements_analyte.keys()):
            value = float(measurements_analyte[sample])
            denominator_values.append(value)
    numerator_mean = statistics.mean(numerator_values)
    denominator_mean = statistics.mean(denominator_values)
    fold = numerator_mean / denominator_mean
    return fold


def calculate_analytes_pairs_folds(
    summary=None,
    group_numerator=None,
    group_denominator=None,
    samples=None,
    measurements=None
):
    """
    Determines the mean fold changes for each analyte between experimental
    groups.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        samples (list<dict<str>>): information about samples from a study
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    # Determine pairs of samples.
    pairs_samples = determine_pairs_samples(
        group_one=group_numerator,
        group_two=group_denominator,
        samples=samples
    )
    # Determine mean fold changes for each analyte.
    summary_novel = []
    for record in summary:
        identifier = record["identifier"]
        fold = calculate_analyte_pairs_fold(
            identifier=identifier,
            group_numerator=group_numerator,
            group_denominator=group_denominator,
            pairs_samples=pairs_samples,
            measurements=measurements
        )
        record["fold"] = fold
        summary_novel.append(record)
    return summary_novel


def determine_pairs_samples(
    group_one=None,
    group_two=None,
    samples=None
):
    """
    Determines the samples for each group for the same patient.

    arguments:
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        samples (list<dict<str>>): information about samples from a study

    raises:

    returns:
        (dict): samples from each group for each pair

    """

    groups = [group_one, group_two]
    pairs_samples = {}
    for sample_cis in samples:
        identifier_cis = sample_cis["identifier"]
        pair_cis = str(sample_cis["pair"])
        group_cis = sample_cis["group"]
        # Determine whether the sample belongs to a relevant group.
        if group_cis in groups:
            # Determine whether the sample belongs to a novel pair.
            if pair_cis not in pairs_samples.keys():
                # Find sample in same pair and other group.
                for sample_trans in samples:
                    identifier_trans = sample_trans["identifier"]
                    pair_trans = str(sample_trans["pair"])
                    group_trans = sample_trans["group"]
                    pair = (pair_trans == pair_cis)
                    group = (
                        (group_trans in groups) and (group_trans != group_cis)
                    )
                    if (pair and group):
                        # Found the other sample for the same patient.
                        break
                pairs_samples[pair_cis] = {
                    group_cis: identifier_cis,
                    group_trans: identifier_trans
                }
    return pairs_samples


def calculate_analyte_pairs_fold(
    identifier=None,
    group_numerator=None,
    group_denominator=None,
    pairs_samples=None,
    measurements=None
):
    """
    Determines the mean fold change between experimental groups for an analyte.

    arguments:
        identifier (str): identifier of an analyte
        group_numerator (str): name of experimental group for numerator
        group_denominator (str): name of experimental group for denominator
        pairs_samples (dict): pairs of samples from both groups
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (float): fold change

    """

    # Find measurements for analyte.
    measurements_analyte = find_analyte_measurements(
        identifier=identifier,
        measurements=measurements
    )
    # Determine fold changes for analyte's measurements.
    folds = []
    for pair in pairs_samples.values():
        sample_numerator = pair[group_numerator]
        sample_denominator = pair[group_denominator]
        if (
            (sample_numerator in measurements_analyte.keys()) and
            (sample_denominator in measurements_analyte.keys())
        ):
            numerator = float(measurements_analyte[sample_numerator])
            denominator = float(measurements_analyte[sample_denominator])
            fold = numerator / denominator
            folds.append(fold)
    mean = statistics.mean(folds)
    return mean


def find_analyte_measurements(identifier=None, measurements=None):
    """
    Finds information about measurements for an analyte.

    arguments:
        identifier (str): identifier of an analyte
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (dict<str>): measurements for an analyte

    """

    def match(record):
        analyte_comparison = utility.convert_string_low_alpha_num(
            record["analyte"]
        )
        identifier_comparison = utility.convert_string_low_alpha_num(
            identifier
        )
        return analyte_comparison == identifier_comparison
    measurements_analyte = utility.find(
        match=match,
        sequence=measurements
    )
    if measurements_analyte is None:
        print("error finding measurements for analyte: " + identifier)
    return measurements_analyte


def calculate_folds_logarithms(records=None):
    """
    Calculates base-2 logarithms of fold changes in records.

    arguments:
        records (list<dict>): information about measurements of analytes

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    records_novel = []
    for record in records:
        fold = record["fold"]
        record["log_fold"] = math.log(fold, 2)
        records_novel.append(record)
    return records_novel


def calculate_analytes_p_values(
    summary=None,
    group_one=None,
    group_two=None,
    samples=None,
    measurements=None
):
    """
    Determines the mean fold changes for each analyte between experimental
    groups.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        samples (list<dict<str>>): information about samples from a study
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    # Determine samples in each group.
    groups_samples = determine_groups_samples(
        group_one=group_numerator,
        group_two=group_denominator,
        samples=samples
    )
    # Determine p-value for each analyte.
    summary_novel = []
    for record in summary:
        identifier = record["identifier"]
        p_value = calculate_analyte_p_value(
            identifier=identifier,
            group_one=group_one,
            group_two=group_two,
            groups_samples=groups_samples,
            measurements=measurements
        )
        record["p_value"] = p_value
        summary_novel.append(record)
    return summary_novel


def calculate_analyte_p_value(
    identifier=None,
    group_one=None,
    group_two=None,
    groups_samples=None,
    measurements=None
):
    """
    Determines the p-value between experimental groups for an analyte.

    arguments:
        identifier (str): identifier of an analyte
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        groups_samples (dict<list<str>>): samples from both groups
        measurements (list<dict<str>>): information measurements from a study

    raises:

    returns:
        (float): p-value

    """

    # Find measurements for analyte.
    measurements_analyte = find_analyte_measurements(
        identifier=identifier,
        measurements=measurements
    )
    # Collect measurements for samples from both groups.
    one_values = []
    for sample in groups_samples[group_one]:
        if (sample in measurements_analyte.keys()):
            value = float(measurements_analyte[sample])
            one_values.append(value)
    two_values = []
    for sample in groups_samples[group_two]:
        if (sample in measurements_analyte.keys()):
            value = float(measurements_analyte[sample])
            two_values.append(value)
    # Determine p-value.
    t_statistic, p_value = scipy.stats.ttest_ind(one_values, two_values)
    return p_value


def calculate_analytes_pairs_p_values(
    summary=None,
    group_one=None,
    group_two=None,
    samples=None,
    measurements=None
):
    """
    Determines the mean fold changes for each analyte between experimental
    groups.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        samples (list<dict<str>>): information about samples from a study
        measurements (list<dict<str>>): information about measurements from a
            study

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    # Determine pairs of samples.
    pairs_samples = determine_pairs_samples(
        group_one=group_one,
        group_two=group_two,
        samples=samples
    )
    # Determine p-value for each analyte.
    summary_novel = []
    for record in summary:
        identifier = record["identifier"]
        p_value = calculate_analyte_pairs_p_value(
            identifier=identifier,
            group_one=group_one,
            group_two=group_two,
            pairs_samples=pairs_samples,
            measurements=measurements
        )
        record["p_value"] = p_value
        summary_novel.append(record)
    return summary_novel


def calculate_analyte_pairs_p_value(
    identifier=None,
    group_one=None,
    group_two=None,
    pairs_samples=None,
    measurements=None
):
    """
    Determines the p-value between experimental groups for an analyte.

    arguments:
        identifier (str): identifier of an analyte
        group_one (str): name of experimental group
        group_two (str): name of experimental group
        pairs_samples (dict): pairs of samples from both groups
        measurements (list<dict<str>>): information measurements from a study

    raises:

    returns:
        (float): p-value

    """

    # Find measurements for analyte.
    measurements_analyte = find_analyte_measurements(
        identifier=identifier,
        measurements=measurements
    )
    # Collect measurements for pairs of samples from both groups.
    groups_measurements = {
        group_one: [],
        group_two: []
    }
    for pair in pairs_samples.values():
        sample_one = pair[group_one]
        sample_two = pair[group_two]
        if (
            (sample_one in measurements_analyte.keys()) and
            (sample_two in measurements_analyte.keys())
        ):
            groups_measurements[group_one].append(
                float(measurements_analyte[sample_one])
            )
            groups_measurements[group_two].append(
                float(measurements_analyte[sample_two])
            )
    # Determine p-values.
    t_statistic, p_value = scipy.stats.ttest_rel(
        groups_measurements[group_one], groups_measurements[group_two]
    )
    return p_value


def filter_analytes_metabolites(
    summary=None
):
    """
    Filter analytes for those that map to metabolites.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes

    raises:

    returns:
        (list<dict<str>>): information about measurements for analytes

    """

    summary_novel = []
    for record in summary:
        metabolites = record["references"]["metabolite"]
        if len(metabolites) > 0:
            summary_novel.append(record)
    return summary_novel


def convert_summary_text(summary=None):
    """
    Converts information about measurements to text format.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes

    returns:
        (list<dict<str>>): information about measurements for analytes

    raises:

    """

    records_text = []
    for record in summary:
        record_text = {
            "identifier": record["identifier"],
            "name": record["name"],
            "reference_hmdb": ";".join(record["references"]["hmdb"]),
            "reference_pubchem": ";".join(record["references"]["pubchem"]),
            "reference_metabolite": (
                ";".join(record["references"]["metabolite"])
            ),
            "fold": record["fold"],
            "log_fold": record["log_fold"],
            "p_value": record["p_value"]
        }
        records_text.append(record_text)
    return records_text


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


def prepare_curation_report(
    summary=None
):
    """
    Prepares a summary report on curation of metabolic sets and entities.

    arguments:
        summary (list<dict<str>>): information about measurements for analytes

    returns:
        (str): report of summary information

    raises:

    """

    # Count measurements.
    count_analytes = len(summary)
    # Count measurements with references to Human Metabolome Database (HMDB).
    count_hmdb = count_records_with_references(
        references=["hmdb"],
        records=summary
    )
    proportion_hmdb = count_hmdb / count_analytes
    percentage_hmdb = round((proportion_hmdb * 100), 2)
    # Count measurements with references to Human Metabolome Database (HMDB).
    count_pubchem = count_records_with_references(
        references=["pubchem"],
        records=summary
    )
    proportion_pubchem = count_pubchem / count_analytes
    percentage_pubchem = round((proportion_pubchem * 100), 2)
    # Count measurements with references to Human Metabolome Database (HMDB).
    count_metab = count_records_with_references(
        references=["metabolite"],
        records=summary
    )
    proportion_metabolite = count_metab / count_analytes
    percent_metabolite = round((proportion_metabolite * 100), 2)
    # Compile information.
    report = textwrap.dedent("""\

        --------------------------------------------------
        curation report

        analytes: {count_analytes}

        measurements with metabolite: {count_metab} ({percent_metabolite} %)
        measurements with PubChem: {count_pubchem} ({percentage_pubchem} %)
        measurements with HMDB: {count_hmdb} ({percentage_hmdb} %)

        --------------------------------------------------
    """).format(
        count_analytes=count_analytes,
        count_hmdb=count_hmdb,
        percentage_hmdb=percentage_hmdb,
        count_pubchem=count_pubchem,
        percentage_pubchem=percentage_pubchem,
        count_metab=count_metab,
        percent_metabolite=percent_metabolite
    )
    # Return information.
    return report


def count_records_with_references(references=None, records=None):
    """
    Counts entities with any of specific references.

    arguments:
        references (list<str>): identifiers of references
        records (list<dict>): information in records

    returns:
        (int): count of records with specific reference

    raises:

    """

    count = 0
    for record in records:
        matches = []
        for reference in references:
            if reference in record["references"].keys():
                if len(record["references"][reference]) > 0:
                    matches.append(True)
        if any(matches):
            count += 1
    return count


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
    path_study_one_pickle = os.path.join(path, "study_one.pickle")
    path_study_one_text = os.path.join(path, "study_one.tsv")
    path_study_two_pickle = os.path.join(path, "study_two.pickle")
    path_study_two_text = os.path.join(path, "study_two.tsv")
    path_study_three_pickle = os.path.join(path, "study_three.pickle")
    path_study_three_text = os.path.join(path, "study_three.tsv")
    path_study_four_pickle = os.path.join(path, "study_four.pickle")
    path_study_four_text = os.path.join(path, "study_four.tsv")
    # Write information to file.
    with open(path_study_one_pickle, "wb") as file_product:
        pickle.dump(information["study_one"]["summary"], file_product)
    utility.write_file_table(
        information=information["study_one"]["summary_text"],
        path_file=path_study_one_text,
        names=information["study_one"]["summary_text"][0].keys(),
        delimiter="\t"
    )
    with open(path_study_two_pickle, "wb") as file_product:
        pickle.dump(information["study_two"]["summary"], file_product)
    utility.write_file_table(
        information=information["study_two"]["summary_text"],
        path_file=path_study_two_text,
        names=information["study_two"]["summary_text"][0].keys(),
        delimiter="\t"
    )
    with open(path_study_three_pickle, "wb") as file_product:
        pickle.dump(information["study_three"]["summary"], file_product)
    utility.write_file_table(
        information=information["study_three"]["summary_text"],
        path_file=path_study_three_text,
        names=information["study_three"]["summary_text"][0].keys(),
        delimiter="\t"
    )
    with open(path_study_four_pickle, "wb") as file_product:
        pickle.dump(information["study_four"]["summary"], file_product)
    utility.write_file_table(
        information=information["study_four"]["summary_text"],
        path_file=path_study_four_text,
        names=information["study_four"]["summary_text"][0].keys(),
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
    # Curate measurements from study zero.
    # Measurements from study one represent metabolites in plasma before and
    # after exercise.
    if False:
        study_zero = curate_measurements_study_zero(
            measurements=source["study_zero"]["measurements"]
        )
    # Curate measurements from study one.
    # Measurements from study two represent metabolites in visceral versus
    # subcutaneous adipose.
    study_one = curate_measurements_study_pairs(
        group_numerator="visceral_fat",
        group_denominator="subcutaneous_fat",
        samples=source["study_one"]["samples"],
        analytes=source["study_one"]["analytes"],
        measurements=source["study_one"]["measurements"],
        hmdb=source["reference"]["hmdb"],
        metabolites=source["reference"]["metabolites"]
    )
    # Curate measurements from study two.
    # Measurements from study two represent metabolites in normal versus tumor
    # lung.
    study_two = curate_measurements_study_pairs(
        group_numerator="tumor",
        group_denominator="normal",
        samples=source["study_two"]["samples"],
        analytes=source["study_two"]["analytes"],
        measurements=source["study_two"]["measurements"],
        hmdb=source["reference"]["hmdb"],
        metabolites=source["reference"]["metabolites"]
    )
    # Curate measurements from study three.
    study_three = curate_measurements_study(
        group_numerator="ischemia",
        group_denominator="normal",
        samples=source["study_three_four"]["samples"],
        analytes=source["study_three_four"]["analytes"],
        measurements=source["study_three_four"]["measurements"],
        hmdb=source["reference"]["hmdb"],
        metabolites=source["reference"]["metabolites"]
    )
    # Curate measurements from study four.
    study_four = curate_measurements_study(
        group_numerator="steatosis",
        group_denominator="normal",
        samples=source["study_three_four"]["samples"],
        analytes=source["study_three_four"]["analytes"],
        measurements=source["study_three_four"]["measurements"],
        hmdb=source["reference"]["hmdb"],
        metabolites=source["reference"]["metabolites"]
    )
    # Compile information.
    information = {
        "study_one": study_one,
        "study_two": study_two,
        "study_three": study_three,
        "study_four": study_four
    }
    #Write product information to file
    write_product(directory=directory, information=information)
