"""
Module to provide template of common structure of modules.

Title:

    template

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
    argparse: Package to interpret user parameters from terminal.
    textwrap: Package to format text.
    csv: Package to organize information in text.
    copy: Package to copy objects.
    pickle: Package to preserve information.
    numpy: Package to calculate with arrays of numbers.
    pandas: Package to organize collections of variables.

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
    University of Utah
    Room 5520C, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of project metabonet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports custom definition of metabolic networks.
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

###############################################################################
# Installation and importation

# Standard

# Relevant

# Custom

#dir()
#importlib.reload()

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
    path_compartments = os.path.join(
        directory, "enhancement_compartments.pickle"
    )
    path_processes = os.path.join(directory, "enhancement_processes.pickle")
    path_reactions = os.path.join(directory, "enhancement_reactions.pickle")
    path_metabolites = os.path.join(
        directory, "enhancement_metabolites.pickle"
    )
    path_curation_compartments = os.path.join(
        directory, "curation_compartments.tsv"
    )
    path_curation_processes = os.path.join(directory, "curation_processes.tsv")
    path_curation_reactions = os.path.join(directory, "curation_reactions.tsv")
    path_curation_metabolites = os.path.join(
        directory, "curation_metabolites.tsv"
    )
    # Read information from file.
    with open(path_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    curation_compartments = utility.read_file_table(
        path_file=path_curation_compartments,
        names=None,
        delimiter="\t"
    )
    curation_processes = utility.read_file_table(
        path_file=path_curation_processes,
        names=None,
        delimiter="\t"
    )
    curation_metabolites = utility.read_file_table(
        path_file=path_curation_metabolites,
        names=None,
        delimiter="\t"
    )
    curation_reactions = utility.read_file_table(
        path_file=path_curation_reactions,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites,
        "curation_compartments": curation_compartments,
        "curation_processes": curation_processes,
        "curation_reactions": curation_reactions,
        "curation_metabolites": curation_metabolites,
    }


def change_compartments(
    changer_compartments=None, compartments_original=None,
    reactions_original=None
):
    """
    Changes information about specific compartments and relevant reactions

    arguments:
        changer_compartments (list<dict<str>>): information to change about
            specific compartments
        compartments_original (dict<dict>): information about processes
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict<dict>>): information about compartments and reactions

    raises:

    """

    # Copy information
    compartments_novel = copy.deepcopy(compartments_original)
    reactions_novel = copy.deepcopy(reactions_original)
    for record in changer_compartments:
        identifier_original = record["identifier_original"]
        identifier_novel = record["identifier_novel"]
        name_original = record["name_original"]
        name_novel = record["name_novel"]
        # Determine method to change information
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel
        if match_identifiers and not match_names:
            # Change name
            if identifier_original in compartments_novel:
                compartments_novel[identifier_original]["name"] = name_novel
        elif match_names and not match_identifiers:
            if identifier_original in compartments_novel:
                # Remove
                del compartments_novel[identifier_original]
                # Remove relevant reactions
                removals = []
                for key, record_reaction in reactions_novel.items():
                    # Determine whether any of reaction's participants are in
                    # the compartment
                    match = determine_reaction_compartment(
                        compartment=identifier_original,
                        reaction=record_reaction
                    )
                    if match:
                        # Remove
                        removals.append(key)
                for removal in removals:
                    del reactions_novel[removal]
    # Compile and return information
    return {
        "compartments": compartments_novel,
        "reactions": reactions_novel
    }


def change_processes(
    changer_processes=None, processes_original=None, reactions_original=None
):
    """
    Changes information about specific processes and relevant reactions

    arguments:
        changer_processes (list<dict<str>>): information to change about
            specific processes
        processes_original (dict<dict>): information about processes
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict<dict>>): information about processes and reactions

    raises:

    """

    # Copy information
    processes_novel = copy.deepcopy(processes_original)
    reactions_novel = copy.deepcopy(reactions_original)
    for record in changer_processes:
        identifier_original = record["identifier_original"]
        identifier_novel = record["identifier_novel"]
        name_original = record["name_original"]
        name_novel = record["name_novel"]
        # Determine method to change information
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel
        if match_identifiers and not match_names:
            # Change name
            if identifier_original in processes_novel:
                processes_novel[identifier_original]["name"] = name_novel
        elif match_names and not match_identifiers:
            if identifier_original in processes_novel:
                # Remove
                del processes_novel[identifier_original]
            if identifier_novel in processes_novel:
                # Replace
                for record_reaction in reactions_novel.values():
                    processes_reaction = record_reaction["processes"]
                    if identifier_original in processes_reaction:
                        for index, process in enumerate(processes_reaction):
                            if process == identifier_original:
                                processes_reaction[index] = identifier_novel
    # Compile and return information
    return {
        "processes": processes_novel,
        "reactions": reactions_novel
    }


def change_metabolites(
    changer_metabolites=None, metabolites_original=None, reactions_original=None
):
    """
    Changes information about specific metabolites and relevant reactions

    arguments:
        changer_metabolites (list<dict<str>>): information to change about
            specific metabolites
        metabolites_original (dict<dict>): information about metabolites
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict<dict>>): information about metabolites and reactions

    raises:

    """

    # Copy information
    metabolites_novel = copy.deepcopy(metabolites_original)
    reactions_novel = copy.deepcopy(reactions_original)
    for record in changer_metabolites:
        identifier_original = record["identifier_original"]
        identifier_novel = record["identifier_novel"]
        name_original = record["name_original"]
        name_novel = record["name_novel"]
        # Determine method to change information
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel
        if match_identifiers and not match_names:
            # Change name
            if identifier_original in metabolites_novel:
                metabolites_novel[identifier_original]["name"] = name_novel
        elif match_names and not match_identifiers:
            if identifier_original in metabolites_novel:
                # Remove
                del metabolites_novel[identifier_original]
            if identifier_novel in metabolites_novel:
                # Replace
                for record_reaction in reactions_novel.values():
                    participants = record_reaction["participants"]
                    for participant in participants:
                        if participant["metabolite"] == identifier_original:
                            participant["metabolite"] = identifier_novel
    # Compile and return information
    return {
        "metabolites": metabolites_novel,
        "reactions": reactions_novel
    }


def change_reactions(changer_reactions=None, reactions_original=None):
    """
    Changes information about specific reactions

    arguments:
        changer_reactions (list<dict<str>>): information to change about
            specific reactions
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    # Copy information
    reactions_novel = copy.deepcopy(reactions_original)
    for record in changer_reactions:
        identifier_original = record["identifier_original"]
        identifier_novel = record["identifier_novel"]
        name_original = record["name_original"]
        name_novel = record["name_novel"]
        # Determine method to change information
        match_identifiers = identifier_original == identifier_novel
        match_names = name_original == name_novel
        if match_identifiers and not match_names:
            # Change name
            if identifier_original in reactions_novel:
                reactions_novel[identifier_original]["name"] = name_novel
        elif match_names and not match_identifiers:
            if identifier_original in reactions_novel:
                # Remove
                del reactions_novel[identifier_original]
    # Return information
    return reactions_novel


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to curate information about metabolic
    entities and sets.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=origin)


    # Change procedures allow custom changes to metabolites and reactions
    # Curate information about compartments.
    compartments = change_compartments(
        changer_compartments=source["changer_compartments"],
        compartments_original=source["compartments"],
        reactions_original=reactions_names
    )
    # Curate information about processes.
    processes = change_processes(
        changer_processes=source["changer_processes"],
        processes_original=source["processes"],
        reactions_original=compartments_change["reactions"]
    )
    # Curate information about reactions.
    reactions = change_reactions(
        changer_reactions=source["changer_reactions"],
        reactions_original=metabolites_change["reactions"]
    )
    # Curate information about metabolites.
    metabolites = change_metabolites(
        changer_metabolites=source["changer_metabolites"],
        metabolites_original=metabolites_references,
        reactions_original=processes_change["reactions"]
    )
