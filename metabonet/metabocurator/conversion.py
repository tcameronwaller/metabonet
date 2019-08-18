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

# The purpose of this procedure is to organize and convert information about
# metabolic sets and entities for further use in the web application.

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
import json

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
    path = os.path.join(directory, "curation")
    path_compartments = os.path.join(path, "compartments.pickle")
    path_processes = os.path.join(path, "processes.pickle")
    path_reactions = os.path.join(path, "reactions.pickle")
    path_metabolites = os.path.join(path, "metabolites.pickle")
    # Read information from file.
    with open(path_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites
    }


def convert_dymetabonet(
    compartments=None, processes=None, reactions=None, metabolites=None
):
    """
    Converts information about metabolic entities and sets to format for web
    applications.

    arguments:
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites

    raises:

    returns:
        (dict<dict<dict>>): information about metabolic entities and sets

    """

    return {
        "compartments": compartments,
        "processes": processes,
        "metabolites": metabolites,
        "reactions": reactions,
    }


def convert_compartments_text(compartments=None):
    """
    Converts information about compartments to text format.

    arguments:
        compartments (dict<dict>): information about compartments

    raises:

    returns:
        (list<dict>): information about compartments

    """

    records = []
    for compartment in compartments.values():
        record = {
            "identifier": compartment["identifier"],
            "name": compartment["name"]
        }
        records.append(record)
    return records


def convert_processes_text(processes=None):
    """
    Converts information about processes to text format.

    arguments:
        processes (dict<dict>): information about processes

    raises:

    returns:
        (list<dict>): information about processes

    """

    records = []
    for process in processes.values():
        record = {
            "identifier": process["identifier"],
            "name": process["name"]
        }
        records.append(record)
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
        # Participants.
        compartments = utility.collect_value_from_records(
            key="compartment", records=reaction["participants"]
        )
        compartments_unique = utility.collect_unique_elements(
            elements_original=compartments
        )
        metabolites = utility.collect_value_from_records(
            key="metabolite", records=reaction["participants"]
        )
        metabolites_unique = utility.collect_unique_elements(
            elements_original=metabolites
        )
        # Transports.
        transport_metabolites = utility.collect_value_from_records(
            key="metabolite", records=reaction["transports"]
        )
        transport_compartments = utility.collect_values_from_records(
            key="compartments", records=reaction["transports"]
        )
        transport_compartments_unique = utility.collect_unique_elements(
            elements_original=transport_compartments
        )
        # Compile information.
        record = {
            "identifier": reaction["identifier"],
            "name": reaction["name"],
            "equation": reaction["equation"],
            "metabolites": ";".join(metabolites_unique),
            "compartments": ";".join(compartments_unique),
            "processes": ";".join(reaction["processes"]),
            "reversibility": reaction["reversibility"],
            "conversion": reaction["conversion"],
            "dispersal": reaction["dispersal"],
            "transport": reaction["transport"],
            "transport_metabolites": ";".join(transport_metabolites),
            "transport_compartments": ";".join(transport_compartments_unique),
            "replication": reaction["replication"],
            "replicates": ";".join(reaction["replicates"]),
            "reference_metanetx": ";".join(reaction["references"]["metanetx"]),
            "reference_recon2m2": ";".join(reaction["references"]["recon2m2"]),
            "reference_gene": ";".join(reaction["references"]["gene"]),
            "reference_enzyme": ";".join(reaction["references"]["enzyme"]),
            "reference_kegg": ";".join(reaction["references"]["kegg"]),
            "reference_reactome": ";".join(reaction["references"]["reactome"]),
            "reference_metacyc": ";".join(reaction["references"]["metacyc"]),
            "reference_bigg": ";".join(reaction["references"]["bigg"]),
            "reference_rhea": ";".join(reaction["references"]["rhea"]),
            "reference_sabiork": ";".join(reaction["references"]["sabiork"]),
            "reference_seed": ";".join(reaction["references"]["seed"])
        }
        records.append(record)
    return records


def convert_reactions_export_text(
    reactions=None,
    metabolites=None,
    compartments=None,
    processes=None,
):
    """
    Converts information about reactions to text format.

    Converts identifiers of metabolites, compartments, and processes to names.

    arguments:
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes


    returns:
        (list<dict>): information about reactions

    raises:

    """

    records = []
    for reaction in reactions.values():
        # Participants.
        # Write a function to compose identifier (name) human readable...
        # Compartments
        compartments_identifiers = utility.collect_value_from_records(
            key="compartment", records=reaction["participants"]
        )
        compartments_identifiers_unique = utility.collect_unique_elements(
            elements_original=compartments_identifiers
        )
        compartments_names = utility.collect_values_from_records_in_reference(
            key="name",
            identifiers=compartments_identifiers_unique,
            reference=compartments,
        )
        # Processes
        processes_names = utility.collect_values_from_records_in_reference(
            key="name",
            identifiers=reaction["processes"],
            reference=processes,
        )
        # Metabolites
        reactants_identifiers = utility.collect_reaction_participants_value(
            key="metabolite",
            criteria={"roles": ["reactant"]},
            participants=reaction["participants"]
        )
        reactants_names = utility.collect_values_from_records_in_reference(
            key="name",
            identifiers=reactants_identifiers,
            reference=metabolites,
        )
        products_identifiers = utility.collect_reaction_participants_value(
            key="metabolite",
            criteria={"roles": ["product"]},
            participants=reaction["participants"]
        )
        products_names = utility.collect_values_from_records_in_reference(
            key="name",
            identifiers=products_identifiers,
            reference=metabolites,
        )
        # Compile information.
        record = {
            "identifier": reaction["identifier"],
            "name": reaction["name"],
            "reactants": "; ".join(reactants_names),
            "products": "; ".join(products_names),
            "compartments": "; ".join(compartments_names),
            "processes": ";".join(processes_names),

            "reversibility": reaction["reversibility"],
            "reference_metanetx": (
                "; ".join(reaction["references"]["metanetx"])
            ),
            "reference_recon2m2": (
                "; ".join(reaction["references"]["recon2m2"])
            ),
            "reference_gene": "; ".join(reaction["references"]["gene"]),
            "reference_enzyme": "; ".join(reaction["references"]["enzyme"]),
            "reference_kegg": "; ".join(reaction["references"]["kegg"]),
            "reference_reactome": (
                "; ".join(reaction["references"]["reactome"])
            ),
            "reference_metacyc": "; ".join(reaction["references"]["metacyc"]),
            "reference_bigg": "; ".join(reaction["references"]["bigg"]),
        }
        records.append(record)
    return records


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
        record = {
            "identifier": metabolite["identifier"],
            "name": metabolite["name"],
            "formula": metabolite["formula"],
            "mass": metabolite["mass"],
            "charge": metabolite["charge"],
            "reference_metanetx":
                ";".join(metabolite["references"]["metanetx"]),
            "reference_hmdb": ";".join(metabolite["references"]["hmdb"]),
            "reference_pubchem": ";".join(metabolite["references"]["pubchem"]),
            "reference_chebi": ";".join(metabolite["references"]["chebi"]),
            "reference_bigg": ";".join(metabolite["references"]["bigg"]),
            "reference_kegg": ";".join(metabolite["references"]["kegg"]),
            "reference_metacyc": ";".join(metabolite["references"]["metacyc"]),
            "reference_reactome":
                ";".join(metabolite["references"]["reactome"]),
            "reference_lipidmaps":
                ";".join(metabolite["references"]["lipidmaps"]),
            "reference_sabiork": ";".join(metabolite["references"]["sabiork"]),
            "reference_seed": ";".join(metabolite["references"]["seed"]),
            "reference_slm": ";".join(metabolite["references"]["slm"]),
            "reference_envipath":
                ";".join(metabolite["references"]["envipath"]),
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
    path = os.path.join(directory, "model")
    utility.confirm_path_directory(path)
    path_dymetabonet = os.path.join(path, "dymetabonet.json")
    path_compartments = os.path.join(path, "compartments.pickle")
    path_processes = os.path.join(path, "processes.pickle")
    path_reactions = os.path.join(path, "reactions.pickle")
    path_metabolites = os.path.join(path, "metabolites.pickle")
    path_compartments_text = os.path.join(path, "compartments.tsv")
    path_processes_text = os.path.join(path, "processes.tsv")
    path_reactions_text = os.path.join(path, "reactions.tsv")
    path_metabolites_text = os.path.join(path, "metabolites.tsv")
    # Write information to file.
    with open(path_dymetabonet, "w") as file_product:
        json.dump(information["dymetabonet"], file_product)
    with open(path_compartments, "wb") as file_product:
        pickle.dump(information["compartments"], file_product)
    with open(path_processes, "wb") as file_product:
        pickle.dump(information["processes"], file_product)
    with open(path_reactions, "wb") as file_product:
        pickle.dump(information["reactions"], file_product)
    with open(path_metabolites, "wb") as file_product:
        pickle.dump(information["metabolites"], file_product)
    utility.write_file_table(
        information=information["compartments_text"],
        path_file=path_compartments_text,
        names=information["compartments_text"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["processes_text"],
        path_file=path_processes_text,
        names=information["processes_text"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["reactions_text"],
        path_file=path_reactions_text,
        names=information["reactions_text"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["metabolites_text"],
        path_file=path_metabolites_text,
        names=information["metabolites_text"][0].keys(),
        delimiter="\t"
    )


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to convert information about metabolic
    entities and sets to versatile formats.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Convert information for export to DyMetaboNet.
    dymetabonet = convert_dymetabonet(
        compartments=source["compartments"],
        processes=source["processes"],
        reactions=source["reactions"],
        metabolites=source["metabolites"]
    )
    # Convert information for export in text.
    compartments_text = convert_compartments_text(
        compartments=source["compartments"]
    )
    processes_text = convert_processes_text(
        processes=source["processes"]
    )
    reactions_text = convert_reactions_text(
        reactions=source["reactions"]
    )
    metabolites_text = convert_metabolites_text(
        metabolites=source["metabolites"]
    )
    # Compile information.
    information = {
        "dymetabonet": dymetabonet,
        "compartments": source["compartments"],
        "processes": source["processes"],
        "metabolites": source["metabolites"],
        "reactions": source["reactions"],
        "compartments_text": compartments_text,
        "processes_text": processes_text,
        "reactions_text": reactions_text,
        "metabolites_text": metabolites_text,
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
