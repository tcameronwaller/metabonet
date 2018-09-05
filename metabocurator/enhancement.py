"""
Enhance information about metabolic sets and entities from MetaNetX.

Title:
    enhancement

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

# Packages and modules from third parties

#import numpy
#import pandas
#import scipy

# Packages and modules from local source

import utility
import extraction

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
    path_hmdb = os.path.join(directory, "hmdb_metabolites_references.pickle")
    path = os.path.join(directory, "extraction")
    path_compartments = os.path.join(path, "extraction_compartments.pickle")
    path_processes = os.path.join(path, "extraction_processes.pickle")
    path_reactions = os.path.join(path, "extraction_reactions.pickle")
    path_metabolites = os.path.join(path, "extraction_metabolites.pickle")
    # Read information from file.
    with open(path_hmdb, "rb") as file_source:
        metabolites_references = pickle.load(file_source)
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
        "metabolites_references": metabolites_references,
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites
    }


def enhance_metabolites(
    metabolites_original=None, metabolites_references=None
):
    """
    Enhances information about metabolites

    arguments:
        metabolites_original (dict<dict>): information about metabolites
        metabolites_references (dict<dict>): metabolites' references from Human
            Metabolome Database (HMDB)

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
            metabolites_references=metabolites_references
        )
        # Compile information
        metabolites_novel[metabolite_novel["identifier"]] = metabolite_novel
    return metabolites_novel


def enhance_metabolite(
    metabolite_original=None, metabolites_references=None
):
    """
    Enhances information about a metabolite

    arguments:
        metabolite_original (dict): information about a metabolite
        metabolites_references (dict<dict>): metabolites' references from Human
            Metabolome Database (HMDB)

    returns:
        (dict): information about a metabolite

    raises:

    """

    # Copy information
    metabolite_novel = copy.deepcopy(metabolite_original)
    # Determine supplemental references from entries in HMDB
    references_supplemental = enhance_metabolite_references(
        references_original=metabolite_novel["references"],
        metabolites_references=metabolites_references
    )
    metabolite_novel["references"] = references_supplemental
    return metabolite_novel


def enhance_metabolite_references(
    references_original=None, metabolites_references=None
):
    """
    Enhances information about a metabolite by including references from HMDB

    arguments:
        references_original (dict): references about a metabolite
        metabolites_references (dict<dict>): metabolites' references from Human
            Metabolome Database (HMDB)

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
        metabolites_references=metabolites_references
    )
    # Extract references from entries in HMDB
    hmdb_references = collect_hmdb_entries_references(
        keys=hmdb_keys, metabolites_references=metabolites_references
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
    identifiers=None, metabolites_references=None
):
    """
    Filters entries from HMDB by their identifiers

    arguments:
        identifiers (list<str>): identifiers by which to find entries in HMDB
        metabolites_references (dict<dict>): metabolites' references from Human
            Metabolome Database (HMDB)

    returns:
        (list<str>): keys of entries in HMDB

    raises:

    """

    keys = []
    for key, record in metabolites_references.items():
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
    keys=None, metabolites_references=None
):
    """
    Extracts references from entries in HMDB

    arguments:
        keys (list<str>): keys of entries in HMDB
        metabolites_references (dict<dict>): metabolites' references from Human
            Metabolome Database (HMDB)

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
        record = metabolites_references[key]
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


def filter_reactions(reactions_original=None):
    """
    Filters reactions by relevance to contextual metabolic network.

    arguments:
        reactions (dict<dict>): information about reactions

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    reactions_novel = {}
    for reaction in reactions_original.values():
        match, record = filter_reaction(reaction=reaction)
        if match:
            reactions_novel[reaction["identifier"]] = record
    return reactions_novel


def filter_reaction(reaction=None):
    """
    Filters a reaction by relevance to contextual metabolic network.

    arguments:
        reaction (dict): information about a reaction

    returns:
        (tuple<bool, dict>): whether reaction passes filters, and report

    raises:

    """

    # Name.
    # Biomass and protein degradation are irrelevant.
    name = (
        reaction["name"] == "Generic human biomass reaction" or
        reaction["name"] == "Protein degradation"
    )
    # Compartment.
    # Boundary is irrelevant.
    compartments = utility.collect_value_from_records(
        key="compartment", records=reaction["participants"]
    )
    compartment = utility.match_string_in_list(
        string="BOUNDARY", list=compartments
    )
    # Metabolite.
    # Biomass is irrelevant.
    metabolites = utility.collect_value_from_records(
        key="metabolite", records=reaction["participants"]
    )
    metabolite = utility.match_string_in_list(
        string="BIOMASS", list=metabolites
    )
    # Reference.
    # MetaNetX reaction MNXR01 is for a meaningless proton exchange.
    reference = reaction["references"]["metanetx"][0] == "MNXR01"
    # Determine whether reaction passes filters.
    filter = name or compartment or metabolite or reference
    # Prepare report.
    record = {
        "identifier_original": reaction["identifier"],
        "identifier_novel": "null",
        "name_original": reaction["name"],
        "name_novel": reaction["name"],
        "custom": False,
        "name": name,
        "compartment": compartment,
        "metabolite": metabolite,
        "reference": reference
    }
    # Return information.
    return (filter, record)


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
    path = os.path.join(directory, "enhancement")
    utility.confirm_path_directory(path)
    path_compartments = os.path.join(path, "enhancement_compartments.pickle")
    path_processes = os.path.join(path, "enhancement_processes.pickle")
    path_reactions = os.path.join(path, "enhancement_reactions.pickle")
    path_metabolites = os.path.join(path, "enhancement_metabolites.pickle")
    path_metabolites_report = os.path.join(path, "metabolites_report.tsv")
    path_reactions_report = os.path.join(path, "reactions_report.tsv")
    path_reactions_filter = os.path.join(path, "reactions_filter.tsv")
    # Write information to file.
    with open(path_compartments, "wb") as file_product:
        pickle.dump(information["compartments"], file_product)
    with open(path_processes, "wb") as file_product:
        pickle.dump(information["processes"], file_product)
    with open(path_reactions, "wb") as file_product:
        pickle.dump(information["reactions"], file_product)
    with open(path_metabolites, "wb") as file_product:
        pickle.dump(information["metabolites"], file_product)
    utility.write_file_table(
        information=information["metabolites_report"],
        path_file=path_metabolites_report,
        names=information["metabolites_report"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["reactions_report"],
        path_file=path_reactions_report,
        names=information["reactions_report"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["reactions_filter"],
        path_file=path_reactions_filter,
        names=information["reactions_filter"][0].keys(),
        delimiter="\t"
    )


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
    source = read_source(directory=directory)
    # Enhance metabolites' references.
    metabolites = enhance_metabolites(
        metabolites_original=source["metabolites"],
        metabolites_references=source["metabolites_references"]
    )
    # Include information about reactions' behavior.
    reactions_behavior = include_reactions_behaviors(
        reactions_original=source["reactions"]
    )
    # Include transport reactions in processes.
    reactions_process = include_reactions_transport_processes(
        reactions_original=reactions_behavior
    )
    # Include information about reactions' replicates.
    reactions_replication = include_reactions_replications(
        reactions_original=reactions_process
    )
    # Prepare reports of information for review.
    metabolites_report = extraction.prepare_report_metabolites(
        metabolites=metabolites
    )
    reactions_report = extraction.prepare_report_reactions(
        reactions=reactions_replication
    )
    # Filter reactions.
    reactions_filter = filter_reactions(
        reactions_original=reactions_replication
    )
    # Compile information.
    information = {
        "compartments": source["compartments"],
        "processes": source["processes"],
        "metabolites": metabolites,
        "reactions": reactions_replication,
        "metabolites_report": metabolites_report,
        "reactions_report": reactions_report,
        "reactions_filter": list(reactions_filter.values())
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
