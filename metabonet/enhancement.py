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
    path_file_hmdb = os.path.join(directory, "hmdb_metabolites.xml")
    path_file_recon2m2 = os.path.join(directory, "recon2m2.xml")
    path_file_compartments = os.path.join(
        directory, "extraction_compartments.pickle"
    )
    path_file_processes = os.path.join(
        directory, "extraction_processes.pickle"
    )
    path_file_metabolites = os.path.join(
        directory, "extraction_metabolites.pickle"
    )
    path_file_reactions = os.path.join(
        directory, "extraction_reactions.pickle"
    )
    path_file_changer_compartments = os.path.join(
        directory, "curation_enhancement_compartments.tsv"
    )
    path_file_changer_processes = os.path.join(
        directory, "curation_enhancement_processes.tsv"
    )
    path_file_changer_metabolites = os.path.join(
        directory, "curation_enhancement_metabolites.tsv"
    )
    path_file_changer_reactions = os.path.join(
        directory, "curation_enhancement_reactions.tsv"
    )
    # Read information from file
    recon2m2 = et.parse(path_file_recon2m2)
    hmdb = et.parse(path_file_hmdb)
    with open(path_file_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_file_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_file_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    with open(path_file_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    changer_compartments = utility.read_file_table(
        path_file=path_file_changer_compartments,
        names=None,
        delimiter="\t"
    )
    changer_processes = utility.read_file_table(
        path_file=path_file_changer_processes,
        names=None,
        delimiter="\t"
    )
    changer_metabolites = utility.read_file_table(
        path_file=path_file_changer_metabolites,
        names=None,
        delimiter="\t"
    )
    changer_reactions = utility.read_file_table(
        path_file=path_file_changer_reactions,
        names=None,
        delimiter="\t"
    )
    # Compile and return information
    return {
        "recon2m2": recon2m2,
        "hmdb": hmdb,
        "compartments": compartments,
        "processes": processes,
        "metabolites": metabolites,
        "reactions": reactions,
        "changer_compartments": changer_compartments,
        "changer_processes": changer_processes,
        "changer_metabolites": changer_metabolites,
        "changer_reactions": changer_reactions
    }


def extract_hmdb_metabolites_references(hmdb=None):
    """
    Extracts metabolites' references from Human Metabolome Database (HMDB)

    arguments:
        hmdb (object): content from HMDB in XML

    returns:
        (dict<dict>): metabolites' references from HMDB

    raises:

    """

    # Interpret content
    reference = utility.interpret_content_hmdb(content=hmdb)
    # Extract information for metabolites
    metabolites_references = {}
    for metabolite in reference["metabolites"].findall(
        "base:metabolite", reference["space"]
    ):
        # HMDB identifiers
        identifier_hmdb_primary = metabolite.find(
            "base:accession", reference["space"]
        ).text
        identifiers_hmdb_secondary = metabolite.find(
            "base:secondary_accessions", reference["space"]
        )
        identifiers_hmdb = [identifier_hmdb_primary]
        for identifier in identifiers_hmdb_secondary.findall(
            "base:accession", reference["space"]
        ):
            identifiers_hmdb.append(identifier.text)
        # Name
        name = metabolite.find("base:name", reference["space"]).text
        # PubChem identifier
        identifier_pubchem = metabolite.find(
            "base:pubchem_compound_id", reference["space"]
        ).text
        # Chemical Entities of Biological Interest (ChEBI) identifier
        identifier_chebi = metabolite.find(
            "base:chebi_id", reference["space"]
        ).text
        # Kyoto Encyclopedia of Genes and Genomes (KEGG) identifier
        identifier_kegg = metabolite.find(
            "base:kegg_id", reference["space"]
        ).text
        # Compile information
        record = {
            "identifier": identifier_hmdb_primary,
            "name": name,
            "hmdb": identifiers_hmdb,
            "pubchem": identifier_pubchem,
            "chebi": identifier_chebi,
            "kegg": identifier_kegg
        }
        metabolites_references[identifier_hmdb_primary] = record
    # Return information
    return metabolites_references


def extract_recon2m2_reactions_names(recon2m2=None):
    """
    Extracts reactions' names from Recon 2M.2

    arguments:
        recon2m2 (object): content from Recon 2M.2 in SBML

    returns:
        (dict<str>): names of reactions from Recon 2M.2

    raises:

    """

    # Copy and interpret content
    reference = utility.copy_interpret_content_recon2m2(content=recon2m2)
    reactions_names = {}
    for reaction in reference["reactions"].findall(
        "version:reaction", reference["space"]
    ):
        identifier_recon2m2 = reaction.attrib["id"]
        name = reaction.attrib["name"]
        reactions_names[identifier_recon2m2] = name
    # Return content with changes
    return reactions_names


def enhance_metabolites(
    metabolites_original=None, hmdb_metabolites_references=None
):
    """
    Enhances information about metabolites

    arguments:
        metabolites_original (dict<dict>): information about metabolites
        hmdb_metabolites_references (dict<dict>): metabolites' references from
            HMDB

    returns:
        (dict<dict>): information about metabolites

    raises:

    """

    metabolites_novel = {}
    # Iterate on records for metabolites
    for key, record in metabolites_original.items():
        # Enhance information about metabolite
        metabolite_novel = enhance_metabolite(
            metabolite_original=record,
            hmdb_metabolites_references=hmdb_metabolites_references
        )
        # Compile information
        metabolites_novel[metabolite_novel["identifier"]] = metabolite_novel
    return metabolites_novel


def enhance_metabolite(
    metabolite_original=None, hmdb_metabolites_references=None
):
    """
    Enhances information about a metabolite

    arguments:
        metabolite_original (dict): information about a metabolite
        hmdb_metabolites_references (dict<dict>): metabolites' references from
            HMDB

    returns:
        (dict): information about a metabolite

    raises:

    """

    # Copy information
    metabolite_novel = copy.deepcopy(metabolite_original)
    # Determine supplemental references from entries in HMDB
    references_supplemental = enhance_metabolite_references(
        references_original=metabolite_novel["references"],
        hmdb_metabolites_references=hmdb_metabolites_references
    )
    metabolite_novel["references"] = references_supplemental
    return metabolite_novel


def enhance_metabolite_references(
    references_original=None, hmdb_metabolites_references=None
):
    """
    Enhances information about a metabolite by including references from HMDB

    arguments:
        references_original (dict): references about a metabolite
        hmdb_metabolites_references (dict<dict>): metabolites' references from
            HMDB

    returns:
        (dict): references about a metabolite

    raises:

    """

    # Copy information
    references_novel = copy.deepcopy(references_original)
    # Each metabolite's record can have references to multiple entries in
    # HMDB
    metabolite_references_hmdb = references_novel["hmdb"]
    # Find entries from HMDB that match
    hmdb_keys = filter_hmdb_entries_identifiers(
        identifiers=metabolite_references_hmdb,
        hmdb_metabolites_references=hmdb_metabolites_references
    )
    # Extract references from entries in HMDB
    hmdb_references = collect_hmdb_entries_references(
        keys=hmdb_keys, hmdb_metabolites_references=hmdb_metabolites_references
    )
    # Combine supplemental references to original references
    references_novel["hmdb"] = utility.collect_unique_elements(
        hmdb_references["hmdb"]
    )
    references_novel["pubchem"] = utility.collect_unique_elements(
        hmdb_references["pubchem"]
    )
    references_novel["chebi"] = utility.collect_unique_elements(
        references_original["chebi"] + hmdb_references["chebi"]
    )
    references_novel["kegg"] = utility.collect_unique_elements(
        references_original["kegg"] + hmdb_references["kegg"]
    )
    return references_novel


def filter_hmdb_entries_identifiers(
    identifiers=None, hmdb_metabolites_references=None
):
    """
    Filters entries from HMDB by their identifiers

    arguments:
        identifiers (list<str>): identifiers by which to find entries in HMDB
        hmdb_metabolites_references (dict<dict>): metabolites' references from
            HMDB

    returns:
        (list<str>): keys of entries in HMDB

    raises:

    """

    keys = []
    for key, record in hmdb_metabolites_references.items():
        hmdb_entry_identifiers = record["hmdb"]
        # Determine whether any of entry's identifiers match the metabolite's
        # references
        checks = []
        for identifier in identifiers:
            check = identifier in hmdb_entry_identifiers
            checks.append(check)
        if any(checks):
            # The entry matches the metabolite's references
            keys.append(key)
    return keys


def collect_hmdb_entries_references(
    keys=None, hmdb_metabolites_references=None
):
    """
    Extracts references from entries in HMDB

    arguments:
        keys (list<str>): keys of entries in HMDB
        hmdb_metabolites_references (dict<dict>): metabolites' references from
            HMDB

    returns:
        (dict<list<str>>): references

    raises:

    """

    hmdb = []
    pubchem = []
    chebi = []
    kegg = []
    for key in keys:
        # Only include primary accession identifiers in collection
        record = hmdb_metabolites_references[key]
        hmdb.append(record["identifier"])
        # Only include valid identifiers in the collection
        if (record["pubchem"] is not None) and (len(record["pubchem"]) > 0):
            pubchem.append(record["pubchem"])
        if (record["chebi"] is not None) and (len(record["chebi"]) > 0):
            chebi.append(record["chebi"])
        if (record["kegg"] is not None) and (len(record["kegg"]) > 0):
            kegg.append(record["kegg"])
    # Compile and return information
    return {
        "hmdb": hmdb,
        "pubchem": pubchem,
        "chebi": chebi,
        "kegg": kegg
    }


def enhance_reactions_names(
    recon2m2_reactions_names=None, reactions_original=None
):
    """
    Enhances reactions by including names from Recon 2M.2

    arguments:
        recon2m2_reactions_names (dict<str>): names of reactions from Recon
        2M.2
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    # Copy information
    reactions_novel = copy.deepcopy(reactions_original)
    for key, record in reactions_novel.items():
        identifier_metanetx = record["identifier"]
        # Reconciliation matches some multiple reactions from Recon 2M.2 to a
        # single reaction in MetaNetX
        identifiers_recon2m2 = record["references"]["recon2m2"]
        identifiers_recon2m2_split = identifiers_recon2m2.split(";")
        identifier_recon2m2 = identifiers_recon2m2_split[0]
        if identifier_recon2m2 in recon2m2_reactions_names:
            name = recon2m2_reactions_names[identifier_recon2m2]
        else:
            name = ""
        reactions_novel[key]["name"] = name
    return reactions_novel


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


def determine_reaction_compartment(compartment=None, reaction=None):
    """
    Determines whether any of reaction's participants are in a compartment

    arguments:
        compartment (str): identifier of a compartment
        reaction (dict): information about a reaction

    returns:
        (bool): whether any of reaction's participants are in the compartment

    raises:

    """

    participants = reaction["participants"]
    for participant in participants:
        if participant["compartment"] == compartment:
            return True
    return False


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


def include_reactions_behaviors(reactions_original=None):
    """
    Includes information about reactions' behavior

    arguments:
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    reactions_novel = {}
    for reaction_original in reactions_original.values():
        reaction_novel = include_reaction_behavior(
            reaction_original=reaction_original
        )
        reactions_novel[reaction_novel["identifier"]] = reaction_novel
    return reactions_novel


def include_reaction_behavior(reaction_original=None):
    """
    Includes information about a reaction's behavior

    arguments:
        reaction_original (dict): information about a reaction

    returns:
        (dict): information about a reaction

    raises:

    """

    # Determine whether reaction involves chemical conversion between reactants
    # and products
    conversion = determine_reaction_conversion(reaction=reaction_original)
    # Determine whether reaction involves participation of metabolites in
    # multiple compartments
    dispersal = determine_reaction_dispersal(reaction=reaction_original)
    # Determine reaction's transports
    transports = collect_reaction_transports(reaction=reaction_original)
    # Determine whether reaction involves transport of metabolites between
    # compartments
    transport = len(transports) > 0
    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel["conversion"] = conversion
    reaction_novel["dispersal"] = dispersal
    reaction_novel["transports"] = transports
    reaction_novel["transport"] = transport
    # Return information
    return reaction_novel


def determine_reaction_conversion(reaction=None):
    """
    Determines whether a reaction involves chemical conversion

    arguments:
        reaction (dict): information about a reaction

    returns:
        (bool): whether the reaction involves chemical conversion

    raises:

    """

    reactant_metabolites = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["reactant"]},
        participants=reaction["participants"]
    )
    product_metabolites = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["product"]},
        participants=reaction["participants"]
    )
    return not utility.compare_lists_by_mutual_inclusion(
        list_one=reactant_metabolites,
        list_two=product_metabolites
    )


def determine_reaction_dispersal(reaction=None):
    """
    Determines whether a reaction involves metabolites in multiple compartments

    arguments:
        reaction (dict): information about a reaction

    returns:
        (bool): whether the reaction involves participation of metabolites in
            multiple compartments

    raises:

    """

    compartments = utility.collect_value_from_records(
        key="compartment", records=reaction["participants"]
    )
    return len(compartments) > 0


def collect_reaction_transports(reaction=None):
    """
    Collects information about a reaction's transports

    arguments:
        reaction (dict): information about a reaction

    returns:
        (list<dict>): information about a reaction's transports

    raises:

    """

    metabolites_reactant = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["reactant"]},
        participants=reaction["participants"]
    )
    metabolites_product = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["product"]},
        participants=reaction["participants"]
    )
    # Collect metabolites that participate as both reactants and products
    metabolites = utility.filter_common_elements(
        list_one=metabolites_product, list_two=metabolites_reactant
    )
    transports = []
    for metabolite in metabolites:
        # Determine metabolite's compartments as reactant and product
        compartments_reactant = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "roles": ["reactant"]
            },
            participants=reaction["participants"]
        )
        compartments_product = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "roles": ["product"]
            },
            participants=reaction["participants"]
        )
        # Determine whether there is a difference between the metabolite's
        # compartments as reactant and product
        transport = not utility.compare_lists_by_mutual_inclusion(
            list_one=compartments_reactant,
            list_two=compartments_product
        )
        if transport:
            compartments = compartments_reactant + compartments_product
            compartments_unique = utility.collect_unique_elements(
                elements_original=compartments
            )
            record = {
                "metabolite": metabolite,
                "compartments": compartments_unique
            }
            transports.append(record)
    return transports


def filter_reaction_participants(criteria=None, participants=None):
    """
    Filters a reaction's participants by multiple criteria

    arguments:
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<dict>): information about a reaction's participants

    raises:

    """

    def match(participant):
        if "metabolites" in criteria:
            match_metabolite = participant["metabolite"] in criteria["metabolites"]
        else:
            match_metabolite = True
        if "compartments" in criteria:
            match_compartment = participant["compartment"] in criteria["compartments"]
        else:
            match_compartment = True
        if "roles" in criteria:
            match_role = participant["role"] in criteria["roles"]
        else:
            match_role = True
        return match_metabolite and match_compartment and match_role
    return list(filter(match, participants))


def collect_reaction_participants_value(
    key=None, criteria=None, participants=None
):
    """
    Collects a value from a reaction's specific participants

    arguments:
        key (str): key of value to collect from each participant
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<str>): values from a reaction's participants

    raises:

    """

    participants_match = filter_reaction_participants(
        criteria=criteria, participants=participants
    )
    return utility.collect_value_from_records(
        key=key, records=participants_match
    )


def include_reactions_transport_processes(reactions_original=None):
    """
    Includes information about reactions' transport processes

    arguments:
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    # Collect information about compartments in which each metabolite
    # participates in each process
    processes_dispersal = collect_processes_metabolites_compartments(
        reactions=reactions_original
    )
    # Filter information for prospective transports of metabolites between
    # compartments in processes
    processes_transports = filter_processes_transports(
        processes_dispersal=processes_dispersal
    )
    reactions_novel = {}
    for reaction_original in reactions_original.values():
        reaction_novel = include_reaction_transport_processes(
            reaction_original=reaction_original,
            processes_transports=processes_transports
        )
        reactions_novel[reaction_novel["identifier"]] = reaction_novel
    return reactions_novel


def collect_processes_metabolites_compartments(reactions=None):
    """
    Collects information about processes' metabolites and compartments

    arguments:
        reactions (dict<dict>): information about reactions

    returns:
        (dict<dict<list<str>>>): information about compartments in which each
            metabolite participates in each process

    raises:

    """

    collection = {}
    for reaction in reactions.values():
        processes = reaction["processes"]
        for process in processes:
            if process not in collection:
                collection[process] = {}
            metabolites = collect_reaction_participants_value(
                key="metabolite",
                criteria={},
                participants=reaction["participants"]
            )
            for metabolite in metabolites:
                if metabolite not in collection[process]:
                    collection[process][metabolite] = []
                compartments = collect_reaction_participants_value(
                    key="compartment",
                    criteria={"metabolites": [metabolite]},
                    participants=reaction["participants"]
                )
                for compartment in compartments:
                    if compartment not in collection[process][metabolite]:
                        collection[process][metabolite].append(compartment)
    return collection


def filter_processes_transports(processes_dispersal=None):
    """
    Collects information about transports in processes

    arguments:
        processes_dispersal (dict<dict<list<str>>>): information about
            compartments in which each metabolite participates in each process

    returns:
        (dict<dict<list<str>>>): information about transports in processes

    raises:

    """

    collection = {}
    for process in processes_dispersal.keys():
        collection[process] = {}
        for metabolite in processes_dispersal[process].keys():
            compartments = processes_dispersal[process][metabolite]
            if len(compartments) > 1:
                collection[process][metabolite] = copy.copy(compartments)
    return collection


def include_reaction_transport_processes(
    reaction_original=None, processes_transports=None
):
    """
    Includes information about a reaction's transport processes

    arguments:
        reaction_original (dict): information about a reaction
        processes_transports (dict<dict<list<str>>>): information about
            transports in processes

    returns:
        (dict): information about a reaction

    raises:

    """

    # Determine processes in which reaction participates by transport
    processes_transport = collect_reaction_transport_processes(
        reaction=reaction_original, processes_transports=processes_transports
    )
    processes_original = reaction_original["processes"]
    processes_total = processes_original + processes_transport
    processes_unique = utility.collect_unique_elements(
        elements_original=processes_total
    )
    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel["processes"] = processes_unique
    # Return information
    return reaction_novel


def collect_reaction_transport_processes(
    reaction=None, processes_transports=None
):
    """
    Collects processes in which a reaction participates by transport

    arguments:
        reaction (dict): information about a reaction
        processes_transports (dict<dict<list<str>>>): information about
            transports in processes

    returns:
        (list<str>): identifiers of processes

    raises:

    """

    transports_reaction = reaction["transports"]
    processes_transport = []
    for process in processes_transports.keys():
        # Determine whether reaction's transports match any of the process'
        # metabolites and compartments
        metabolites_process = processes_transports[process]
        for transport_reaction in transports_reaction:
            metabolite_reaction = transport_reaction["metabolite"]
            compartments_reaction = transport_reaction["compartments"]
            if metabolite_reaction in metabolites_process.keys():
                compartments_process = metabolites_process[metabolite_reaction]
                # Determine whether multiple compartments match between the
                # reaction and the process
                compartments = utility.filter_common_elements(
                    list_one=compartments_reaction,
                    list_two=compartments_process
                )
                if len(compartments) > 1:
                    # Reaction participates in the process by transport
                    processes_transport.append(process)
    return processes_transport


def include_reactions_replications(reactions_original=None):
    """
    Includes information about reactions' replications

    arguments:
        reactions_original (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    # Collect information about replicate reactions
    reactions_replicates = collect_reactions_replicates(
        reactions=reactions_original
    )
    reactions_novel = {}
    for reaction_original in reactions_original.values():
        reaction_novel = include_reaction_replication(
            reaction_original=reaction_original,
            reactions_replicates=reactions_replicates
        )
        reactions_novel[reaction_novel["identifier"]] = reaction_novel
    return reactions_novel


def collect_reactions_replicates(reactions=None):
    """
    Collects information about reactions' replications

    arguments:
        reactions (dict<dict>): information about reactions

    returns:
        (list<dict>): information about reactions' replications

    raises:

    """

    reactions_replicates = []
    for reaction in reactions.values():
        identifier = reaction["identifier"]
        # Collect identifiers of metabolites that participate as reactants and
        # products in the reaction
        reactants = collect_reaction_participants_value(
            key="metabolite",
            criteria={"roles": ["reactant"]},
            participants=reaction["participants"]
        )
        products = collect_reaction_participants_value(
            key="metabolite",
            criteria={"roles": ["product"]},
            participants=reaction["participants"]
        )
        # Determine whether collection includes a record for an identical
        # combination of reactants and products
        index = find_index_reactions_replicates_reactants_products(
            reactants=reactants,
            products=products,
            reactions_replicates=reactions_replicates
        )
        if index == -1:
            # Record does not exist
            # Create novel record
            record = {
                "reactions": [reaction["identifier"]],
                "reactants": reactants,
                "products": products
            }
            reactions_replicates.append(record)
        else:
            # Record exists
            # Include reaction in record
            reactions_replicates[index]["reactions"].append(identifier)
    return reactions_replicates


def find_index_reactions_replicates_reactants_products(
    reactants=None, products=None, reactions_replicates=None
):
    """
    Finds index of a record for replicate reactions by reactants and products

    arguments:
        reactants (list<str>): identifiers of metabolites that participate in a
            reaction as reactants
        products (list<str>): identifiers of metabolites that participate in a
            reaction as products
        reactions_replicates (list<dict>): information about reactions'
            replications

    returns:
        (int): index of record if it exists or -1 if it does not exist

    raises:

    """

    def match (record=None):
        match_reactants = utility.compare_lists_by_mutual_inclusion(
            list_one=reactants,
            list_two=record["reactants"]
        )
        match_products = utility.compare_lists_by_mutual_inclusion(
            list_one=products,
            list_two=record["products"]
        )
        return match_reactants and match_products
    return utility.find_index(match, reactions_replicates)


def include_reaction_replication(
    reaction_original=None, reactions_replicates=None
):
    """
    Includes information about a reaction's replication

    arguments:
        reaction_original (dict): information about a reaction
        reactions_replicates (list<dict>): information about reactions'
            replications

    returns:
        (dict): information about a reaction

    raises:

    """

    # Determine replicate reactions
    index = find_index_reactions_replicates_identifier(
        identifier=reaction_original["identifier"],
        reactions_replicates=reactions_replicates
    )
    if index == -1:
        replicates = []
    else:
        replicates = reactions_replicates[index]["reactions"]
    # Compile information
    reaction_novel = copy.deepcopy(reaction_original)
    reaction_novel["replicates"] = replicates
    reaction_novel["replication"] = len(replicates) > 1
    # Return information
    return reaction_novel


def find_index_reactions_replicates_identifier(
    identifier=None, reactions_replicates=None
):
    """
    Finds index of a record for replicate reactions by identifier

    arguments:
        identifier (str): identifier of a reaction
        reactions_replicates (list<dict>): information about reactions'
            replications

    returns:
        (int): index of record if it exists or -1 if it does not exist

    raises:

    """

    def match (record=None):
        return identifier in record["reactions"]
    return utility.find_index(match, reactions_replicates)


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
    #path_file_hmdb = os.path.join(
    #    directory, "enhancement_hmdb_metabolites_references.csv"
    #)
    #path_file_recon2m2 = os.path.join(
    #    directory, "enhancement_recon2m2_reactions_names.csv"
    #)
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
    # Write information to file
    #utility.write_file_table(
    #    information=information["hmdb_metabolites_references"],
    #    path_file=path_file_hmdb,
    #    names=["identifier", "name", "hmdb", "pubchem", "chebi", "kegg"],
    #    delimiter="\t"
    #)
    #utility.write_file_table(
    #    information=information["recon2m2_reactions_names"],
    #    path_file=path_file_recon2m2,
    #    names=["identifier", "name"],
    #    delimiter="\t"
    #)
    with open(path_file_compartments, "wb") as file_product:
        pickle.dump(information["compartments"], file_product)
    with open(path_file_processes, "wb") as file_product:
        pickle.dump(information["processes"], file_product)
    with open(path_file_metabolites, "wb") as file_product:
        pickle.dump(information["metabolites"], file_product)
    with open(path_file_reactions, "wb") as file_product:
        pickle.dump(information["reactions"], file_product)


###############################################################################
# Procedure


def main():
    """
    This function defines the main activity of the module.
    """

    # Read source information from file
    source = read_source()

    # Extract metabolites' references from Human Metabolome Database
    hmdb_metabolites_references = extract_hmdb_metabolites_references(
        hmdb=source["hmdb"]
    )
    # Extract reactions' names from Recon 2M.2
    recon2m2_reactions_names = extract_recon2m2_reactions_names(
        recon2m2=source["recon2m2"]
    )
    # Enhance metabolites' references
    metabolites_references = enhance_metabolites(
        metabolites_original=source["metabolites"],
        hmdb_metabolites_references=hmdb_metabolites_references
    )
    # Enhance reactions' names
    reactions_names = enhance_reactions_names(
        recon2m2_reactions_names=recon2m2_reactions_names,
        reactions_original=source["reactions"]
    )

    # TODO: Move change procedures to before the enhancement from HMDB...


    # Change procedures allow custom changes to metabolites and reactions
    # Change compartments' information
    compartments_change = change_compartments(
        changer_compartments=source["changer_compartments"],
        compartments_original=source["compartments"],
        reactions_original=reactions_names
    )
    # Change processes information
    processes_change = change_processes(
        changer_processes=source["changer_processes"],
        processes_original=source["processes"],
        reactions_original=compartments_change["reactions"]
    )
    # Change metabolites' information
    metabolites_change = change_metabolites(
        changer_metabolites=source["changer_metabolites"],
        metabolites_original=metabolites_references,
        reactions_original=processes_change["reactions"]
    )
    # Change reactions' information
    reactions_change = change_reactions(
        changer_reactions=source["changer_reactions"],
        reactions_original=metabolites_change["reactions"]
    )



    # Include information about reactions' behavior
    reactions_behavior = include_reactions_behaviors(
        reactions_original=reactions_change
    )
    # Include transport reactions in processes
    reactions_process = include_reactions_transport_processes(
        reactions_original=reactions_behavior
    )
    # Include information about reactions' replicates
    reactions_replication = include_reactions_replications(
        reactions_original=reactions_process
    )
    #Write product information to file
    information = {
        "hmdb_metabolites_references": list(
            hmdb_metabolites_references.values()
        ),
        "recon2m2_reactions_names": list(recon2m2_reactions_names.values()),
        "compartments": compartments_change["compartments"],
        "processes": processes_change["processes"],
        "metabolites": metabolites_change["metabolites"],
        "reactions": reactions_replication
    }
    write_product(information=information)


if __name__ == "__main__":
    main()
