"""
Curate metabolic model Recon 2M.2 for reconciliation to MetaNetX.

Title:
    experiment_group.py

Imports:
    os: This module is from The Python Standard Library. It contains
        difinitions of tools to interact with the operating system.
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

# The purpose of this procedure is to improve reconciliation of information
# about metabolic sets and entities from the Recon 2M.2 model of human
# metabolism to MetaNetX's name space.
# In the process of reconciliation, MetaNetX also exerts some checks and
# changes for quality.

###############################################################################
# Installation and importation of packages and modules


# Packages and modules from the python standard library

import os
#import sys
import shutil
#import importlib
import re
import csv
import xml.etree.ElementTree as et
import copy

# Packages and modules from third parties

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

    # Specify directories and files
    path_model = os.path.join(directory, "recon2m2.xml")
    path_metabolites = os.path.join(
        directory, "reconciliation_metabolites.tsv"
    )
    #path_reactions = os.path.join(
    #    directory, "reconciliation_reactions.tsv"
    #)
    # Read information from file
    model = et.parse(path_model)
    curation_metabolites = utility.read_file_table(
        path_file=path_metabolites,
        names=None,
        delimiter="\t"
    )
    # Compile and return information
    return {
        "model": model,
        "curation_metabolites": curation_metabolites
    }


def remove_model_identifier_prefix(content=None):
    """
    Removes unnecessary prefixes from identifiers for model's entities

    This function removes unnecessary prefixes from identifiers for
    metabolites.

    arguments:
        content (object): content from file in Systems Biology Markup Language
            (XML)

    returns:
        content with changes to identifiers of metabolites

    raises:

    """

    # Copy and interpret content
    reference = utility.copy_interpret_content_recon2m2(content=content)
    # Remove prefixes from identifiers for metabolites
    for metabolite in reference["metabolites"].findall(
    "version:species", reference["space"]
    ):
        # Remove prefix from metabolite's identifier
        novel_identifier = re.sub(r"^M_", "", metabolite.attrib["id"])
        metabolite.attrib["id"] = novel_identifier
        # Search metabolite's annotation
        for description in metabolite.iter(
        "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}description"
        ):
            # Remove prefix from metabolite's identifier
            novel_identifier = re.sub(
                r"^#M_",
                "#",
                description.attrib[
                "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
                ]
            )
            description.attrib[
            "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
            ] = novel_identifier
    # Remove prefixes from identifiers for reactions' metabolites
    for reaction in reference["reactions"].findall(
    "version:reaction", reference["space"]
    ):
        # Search reaction's metabolites
        for metabolite in reaction.iter(
        "{http://www.sbml.org/sbml/level2/version4}speciesReference"
        ):
            # Remove prefix from metabolite's identifier
            novel_identifier = re.sub(r"^M_", "", metabolite.attrib["species"])
            metabolite.attrib["species"] = novel_identifier
    # Return content with changes
    return reference["content"]


def change_model_boundary(content=None):
    """
    Changes annotations for a model's boundary

    This function changes annotations of a model's boundary in compartments,
    metabolites, and reactions.

    arguments:
        content (object): content from file in Systems Biology Markup Language
            (XML)

    returns:
        content with changes to attributes of compartments, metabolites, and
            reactions

    raises:

    """

    # Copy and interpret content
    reference = utility.copy_interpret_content_recon2m2(content=content)
    # Correct compartment for model's boundary
    for compartment in reference["compartments"].findall(
    "version:compartment", reference["space"]
    ):
        # Determine whether compartment is for model's boundary
        if compartment.attrib["id"] == "b":
            compartment.attrib["name"] = "model boundary"
    # Correct metabolites for model's boundary
    for metabolite in reference["metabolites"].findall(
    "version:species", reference["space"]
    ):
        # Determine whether metabolite's compartment is model's boundary
        if "boundary" in metabolite.attrib["id"]:
            novel_identifier = re.sub(
            r"_[eciglmnrx]_boundary", "_b", metabolite.attrib["id"]
            )
            metabolite.attrib["id"] = novel_identifier
            novel_compartment = "b"
            metabolite.attrib["compartment"] = novel_compartment
    # Correct reactions for model's boundary
    for reaction in reference["reactions"].findall(
    "version:reaction", reference["space"]
    ):
        # Search reaction's metabolites
        for metabolite in reaction.iter(
        "{http://www.sbml.org/sbml/level2/version4}speciesReference"
        ):
            # Determine whether metabolite's compartment is model's boundary
            if "boundary" in metabolite.attrib["species"]:
                novel_identifier = re.sub(
                r"_[eciglmnrx]_boundary", "_b", metabolite.attrib["species"]
                )
                metabolite.attrib["species"] = novel_identifier
    # Return content with changes
    return reference["content"]


def change_model_compartment(content=None):
    """
    Changes annotations for a model's compartments

    This function changes annotations of a model's compartments.

    arguments:
        content (object): content from file in Systems Biology Markup Language
            (XML)

    returns:
        content with changes to attributes of compartments, metabolites, and
            reactions

    raises:

    """

    # Copy and interpret content
    reference = utility.copy_interpret_content_recon2m2(content=content)
    # Correct compartment for model's boundary
    for compartment in reference["compartments"].findall(
    "version:compartment", reference["space"]
    ):
        # Determine whether compartment is for model's boundary
        if compartment.attrib["id"] == "e":
            compartment.attrib["name"] = "extracellular region"
    # Return content with changes
    return reference["content"]


def change_model_metabolites_identifiers(
        metabolites_identifiers=None, content=None
):
    """
    Changes metabolites' identifiers

    This function changes metabolites' identifiers according to information
    about translation.

    arguments:
        metabolites_identifiers (list<dict<str>>): translations of metabolites'
            identifiers
        content (object): content from file in Systems Biology Markup Language
            (XML)

    returns:
        content with changes to metabolites' identifiers

    raises:

    """

    # Copy and interpret content
    reference = utility.copy_interpret_content_recon2m2(content=content)
    # Change content for each combination of original and novel identifiers
    for row in metabolites_identifiers:
        # Construct targets to recognize original and novel identifiers
        original_elements = [row["identifier_original"], "_"]
        original_target = "".join(original_elements)
        novel_elements = [row["identifier_novel"], "_"]
        novel_target = "".join(novel_elements)
        # Change identifiers of metabolites
        for metabolite in reference["metabolites"].findall(
        "version:species", reference["space"]
        ):
            # Determine whether to change metabolite's identifier
            if original_target in metabolite.attrib["id"]:
                metabolite.attrib["id"] = metabolite.attrib["id"].replace(
                    original_target, novel_target
                )
        # Change identifiers of reactions' metabolites
        for reaction in reference["reactions"].findall(
        "version:reaction", reference["space"]
        ):
            # Search reaction's metabolites
            for metabolite in reaction.iter(
            "{http://www.sbml.org/sbml/level2/version4}speciesReference"
            ):
                # Determine whether to change metabolite's identifier
                if original_target in metabolite.attrib["species"]:
                    metabolite.attrib["species"] = (
                        metabolite.attrib["species"].replace(
                            original_target, novel_target
                        )
                    )
    # Return content with changes
    return reference["content"]


def write_product(directory=None, content=None):
    """
    Writes product content to file

    arguments:
        directory (str): directory for product files
        content (object): content from file in Systems Biology Markup Language
            (XML)

    returns:

    raises:

    """

    path_file = os.path.join(
        directory, "recon2m2_metanetx_reconciliation.xml"
    )
    # Write information to file
    content.write(path_file, xml_declaration=False)


###############################################################################
# Procedure


def execute_procedure(origin=None, destination=None, clean=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to reconcile the metabolic model to be
    compatible with MetaNetX.

    arguments:
        origin (str): directory of source files
        destination (str): directory for product files
        clean (bool): whether to remove intermediate files

    returns:

    raises:

    """

    print("Executing reconciliation procedure...")

    # Read source information from file
    source = read_source(directory=origin)

    # TODO: should I do the filtering before or after MetaNetX???
    # TODO: before might facilitate or improve reconciliation to MetaNetX... although it's sort of less convenient.
    # TODO: decision... reconciliation should ONLY focus on reconciliation for optimal import to MetaNetX
    # TODO: deal with filtering afterwards.

    # TODO: remove all boundary exchange/demand reactions
    # TODO: remove all boundary metabolites
    # TODO: remove all BIOMASS metabolites
    # TODO: remove all BIOMASS reactions


    if False:
        # Correct content
        model_compartment = change_model_compartment(content=source["model"])
        model_boundary = change_model_boundary(content=content_compartment)
        model_prefix = remove_model_identifier_prefix(content=content_boundary)
        model_identifier = change_model_metabolites_identifiers(
            metabolites_identifiers=source["curation_metabolites"],
            content=model_prefix
        )
        #Write product content to file
        write_product(content=model_identifier)

    pass
