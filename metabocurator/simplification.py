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
import xml.etree.ElementTree as et

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
    path_hmdb = os.path.join(directory, "hmdb_metabolites.xml")
    # Read information from file.
    #hmdb = et.parse(path_hmdb)
    # Compile and return information.
    return {
        "hmdb": path_hmdb
    }


def construct_tag(tag=None, space=None, spaces=None):
    return "{" + spaces[space] + "}" + tag


def extract_hmdb_metabolites_references(hmdb=None):
    """
    Extracts metabolites' references from Human Metabolome Database (HMDB)

    arguments:
        hmdb (str): path to file of HMDB in XML

    raises:

    returns:
        (dict<dict>): metabolites' references from HMDB

    """

    # Collect references to name space.
    spaces = {}
    # Collect references for metabolites.
    metabolites_references = {}
    # Count records.
    count = 0
    for event, element in et.iterparse(
        hmdb, events=('start', 'end', 'start-ns', 'end-ns')
    ):
        if event == "start-ns":
            space = element[0]
            reference = element[1]
            spaces[space] = reference
        if event == "end":
            if element.tag == construct_tag(
                tag="metabolite", space=space, spaces=spaces
            ):
                # Parse complete for a new metabolite.
                # Count.
                count = count + 1
                # HMDB identifiers.
                tag_accession = construct_tag(
                    tag="accession", space=space, spaces=spaces
                )
                hmdb_primary = metabolite.find(tag_accession).text
                tag_secondary = construct_tag(
                    tag="secondary_accessions", space=space, spaces=spaces
                )
                hmdb_secondary = metabolite.find(tag_secondary)
                identifiers_hmdb = [hmdb_primary]
                for identifier in hmdb_secondary.findall(tag_accession):
                    identifiers_hmdb.append(identifier.text)
                # Name.
                #name = element.find("{http://www.hmdb.ca}name").text
                tag_name = construct_tag(
                    tag="name", space=space, spaces=spaces
                )
                name = element.find(tag_name).text
                # PubChem identifier.
                tag_pubchem = construct_tag(
                    tag="pubchem_compound_id", space=space, spaces=spaces
                )
                identifier_pubchem = metabolite.find(tag_pubchem).text
                # Chemical Entities of Biological Interest (ChEBI) identifier.
                tag_chebi = construct_tag(
                    tag="chebi_id", space=space, spaces=spaces
                )
                identifier_chebi = metabolite.find(tag_chebi).text
                # Kyoto Encyclopedia of Genes and Genomes (KEGG) identifier.
                tag_kegg = construct_tag(
                    tag="kegg_id", space=space, spaces=spaces
                )
                identifier_kegg = metabolite.find(tag_kegg).text
                # Compile information
                record = {
                    "identifier": identifier_hmdb_primary,
                    "name": name,
                    "hmdb": identifiers_hmdb,
                    "pubchem": identifier_pubchem,
                    "chebi": identifier_chebi,
                    "kegg": identifier_kegg
                }
                metabolites_references[hmdb_primary] = record
                # Clear memory.
                element.clear()
    # Report.
    print("Extraction complete for " + str(metabolites) + " metabolites!")
    # Return information.
    return metabolites_references


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
    path_pickle = os.path.join(directory, "hmdb_metabolites_references.pickle")
    path_text = os.path.join(directory, "hmdb_metabolites_references.tsv")
    # Write information to file.
    with open(path_pickle, "wb") as file_product:
        pickle.dump(information["references_object"], file_product)
    utility.write_file_table(
        information=information["references_list"],
        path_file=path_text,
        names=["identifier", "name", "hmdb", "pubchem", "chebi", "kegg"],
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

    print("beginning simplification procedure")

    # Read source information from file.
    source = read_source(directory=directory)
    # Extract metabolites' references from Human Metabolome Database.
    metabolites_references = extract_hmdb_metabolites_references(
        hmdb=source["hmdb"]
    )
    if False:
        # Compile information.
        information = {
            "references_object": metabolites_references,
            "references_list": list(metabolites_references.values())
        }
        #Write product information to file
        write_product(directory=directory, information=information)
