"""
Extract information about metabolic sets and entities from MetaNetX.

Title:
    experiment_group.py

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


def read_source():
    """
    Reads and organizes source information from file

    arguments:

    returns:
        (object): source information

    raises:

    """

    # Specify directories and files
    directory = os.path.join(
        os.sep, "media", "tcameronwaller", "primary", "data", "local", "work",
        "project_metabolism", "metabolism_models", "homo_sapiens",
        "recon_2-m-2"
    )
    path_file_compartments = os.path.join(
        directory, "enhancement_compartments.pickle"
    )
    path_file_processes = os.path.join(
        directory, "enhancement_processes.pickle"
    )
    path_file_metabolites = os.path.join(
        directory, "enhancement_metabolites.pickle"
    )
    path_file_reactions = os.path.join(
        directory, "enhancement_reactions.pickle"
    )
    # Read information from file
    with open(path_file_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_file_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_file_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    with open(path_file_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    # Compile and return information
    return {
        "compartments": compartments,
        "processes": processes,
        "metabolites": metabolites,
        "reactions": reactions
    }


def write_product(information=None):
    """
    Writes product information to file

    arguments:
        information (dict): product information

    returns:

    raises:

    """

    # Specify directories and files
    directory = os.path.join(
        os.sep, "media", "tcameronwaller", "primary", "data", "local", "work",
        "project_metabolism", "metabolism_models", "homo_sapiens",
        "recon_2-m-2"
    )
    path_file_metabolism = os.path.join(
        directory, "metabolism_sets_entities_recon2m2.json"
    )
    path_file_metabolites_report = os.path.join(
        directory, "report_metabolites.tsv"
    )
    path_file_reactions_report = os.path.join(
        directory, "report_reactions.tsv"
    )
    # Write information to file
    with open(path_file_metabolism, "w") as file_product:
        json.dump(information["metabolism"], file_product)


###############################################################################
# Procedure


def main():
    """
    This function defines the main activity of the module.
    """

    # Read source information from file
    source = read_source()
    #Write product information to file
    metabolism = {
        "compartments": source["compartments"],
        "processes": source["processes"],
        "metabolites": source["metabolites"],
        "reactions": source["reactions"],
    }
    information = {
        "metabolism": metabolism,
        "metabolites_report": metabolites_report,
        "reactions_report": reactions_report
    }
    write_product(information=information)


if __name__ == "__main__":
    main()
