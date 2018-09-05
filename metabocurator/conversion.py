"""
Convert information about metabolic sets and entities to versatile formats.

Title:
    conversion

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
    Reads and organizes source information from file

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
    path = os.path.join(directory, "conversion")
    utility.confirm_path_directory(path)
    path_dymetabonet = os.path.join(path, "metabolism.json")
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
    # TODO:
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
