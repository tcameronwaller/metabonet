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
import xml.etree.ElementTree as et

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
    path_hmdb = os.path.join(directory, "hmdb_metabolites.xml")
    path = os.path.join(directory, "provision")
    path_midas = os.path.join(path, "midas_original.tsv")
    # Read information from file.
    #hmdb = et.parse(path_hmdb)
    midas_original = utility.read_file_table(
        path_file=path_midas,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "hmdb": path_hmdb,
        "midas_original": midas_original
    }


def construct_tag(tag=None, space=None, spaces=None):
    """
    Constructs complete tag for name space in Extensible Markup Language (XML).

    arguments:
        tag (str): name of element's tag
        space (str): name name space within XML document
        spaces (dict<str>): name spaces within XML document

    raises:

    returns:
        (str): complete name of tag

    """

    return "{" + spaces[space] + "}" + tag


def extract_hmdb_summary(hmdb=None):
    """
    Extracts metabolites' references from Human Metabolome Database (HMDB).

    arguments:
        hmdb (str): path to file of HMDB in XML

    raises:

    returns:
        (dict<dict>): metabolites' references from HMDB

    """

    # Collect references to name space.
    spaces = {}
    # Collect references for metabolites.
    summary_hmdb = {}
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
                # Extract information from record.
                record = extract_hmdb_record_summary(
                    element=element,
                    space=space,
                    spaces=spaces
                )
                summary_hmdb[hmdb_primary] = record
                # Clear memory.
                element.clear()
    # Report.
    print("Extraction complete for " + str(count) + " metabolites!")
    # Return information.
    return summary_hmdb


def extract_hmdb_record_summar(element=None, space=None, spaces=None):

    # HMDB identifiers.
    tag_accession = construct_tag(
        tag="accession", space=space, spaces=spaces
    )
    hmdb_primary = element.find(tag_accession).text
    tag_secondary = construct_tag(
        tag="secondary_accessions", space=space, spaces=spaces
    )
    hmdb_secondary = element.find(tag_secondary)
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
    identifier_pubchem = element.find(tag_pubchem).text
    # Chemical Entities of Biological Interest (ChEBI) identifier.
    tag_chebi = construct_tag(
        tag="chebi_id", space=space, spaces=spaces
    )
    identifier_chebi = element.find(tag_chebi).text
    # Kyoto Encyclopedia of Genes and Genomes (KEGG) identifier.
    tag_kegg = construct_tag(
        tag="kegg_id", space=space, spaces=spaces
    )
    identifier_kegg = element.find(tag_kegg).text
    # Compile and return information.
    record = {
        "identifier": hmdb_primary,
        "name": name,
        "reference_hmdb": ";".join(identifiers_hmdb),
        "reference_pubchem": identifier_pubchem,
        "reference_chebi": identifier_chebi,
        "reference_kegg": identifier_kegg
    }
    return record


def transfer_summary_midas(summary_hmdb=None, midas_original=None):
    """
    Transfers information from Human Metabolome Database (HMDB) to MIDAS
    library.

    arguments:
        summary_hmdb (dict<dict>): information from HMDB
        midas_original (list<dict<str>>): information in MIDAS

    raises:

    returns:
        (list<dict<str>>): information in MIDAS

    """

    midas_novel = []
    for record_midas in midas_original:
        # Interpretation.
        identifier_midas = record_midas["identifier_midas"]
        name_original = record_midas["name"]
        reference_hmdb_original = record_midas["reference_hmdb"]
        reference_kegg_original = record_midas["reference_kegg"]
        # Determine whether MIDAS record matches a HMDB record.
        if len(reference_hmdb_original) > 0:
            # Match MIDAS record to HMDB record.
            keys_hmdb = enhancement.filter_hmdb_entries_identifiers(
                identifiers=[reference_hmdb_original],
                metabolites_references=summary_hmdb
            )
            if len(keys_hmdb) > 1:
                print("found multiple hmdb keys!")
            key_hmdb = keys_hmdb[0]
            record_hmdb = summary_hmdb[key_hmdb]
            record_novel = copy.deepcopy(record_hmdb)
            record_novel["identifier_midas"] = identifier_midas
            record_novel["name_original"] = name_original
            record_novel["reference_kegg_original"] = reference_kegg_original
            midas_novel.append(record_novel)
        else:
            keys_hmdb = list(summary_hmdb.values())[0].keys()
            record_hmdb = {}
            for key in keys_hmdb:
                record_hmdb[key] = "null"
            record_novel = record_hmdb
            record_novel["identifier_midas"] = identifier_midas
            record_novel["name_original"] = name_original
            record_novel["reference_kegg_original"] = reference_kegg_original
            midas_novel.append(record_novel)
    return midas_novel


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
    source = read_source(directory=directory)
    # Extract metabolites' references from Human Metabolome Database.
    summary_hmdb = extract_hmdb_summary(hmdb=source["hmdb"])
    # Transfer information to MIDAS library.
    midas_novel = transfer_summary_midas(
        summary_hmdb=summary_hmdb, midas_original=source["midas_original"]
    )
    # Compile information.
    information = {
        "summary_object": summary_hmdb,
        "summary_list": list(summary_hmdb.values()),
        "midas_novel": midas_novel
    }
    #Write product information to file
    write_product(directory=directory, information=information)
