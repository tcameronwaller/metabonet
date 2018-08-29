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

Author:
    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 5520C, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

"""


###############################################################################
# Notes

# The purpose of this procedure is to extract information about metabolic sets
# and entities from information that MetaNetX exports from its reconciliation.
# This extraction converts information into a convenient format for subsequent
# use and also begins to derive supplemental information.

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
    path_file_genes = os.path.join(
        directory, "recon2m2_metanetx_genes.tsv"
    )
    path_file_compartments = os.path.join(
        directory, "recon2m2_metanetx_compartments.tsv"
    )
    path_file_metabolites = os.path.join(
        directory, "recon2m2_metanetx_metabolites.tsv"
    )
    path_file_reactions = os.path.join(
        directory, "recon2m2_metanetx_reactions.tsv"
    )
    # Read information from file
    #with open(in_file_path_model, "r") as in_file:
    #    content = in_file.read()
    genes = utility.read_file_table(
        path_file=path_file_genes,
        names=["reaction", "genes", "low_bound", "up_bound", "direction"],
        delimiter="\t"
    )
    compartments = utility.read_file_table(
        path_file=path_file_compartments,
        names=["identifier", "name", "source"],
        delimiter="\t"
    )
    metabolites = utility.read_file_table(
        path_file=path_file_metabolites,
        names=[
            "identifier", "name", "source", "formula", "mass", "charge",
            "references"
        ],
        delimiter="\t"
    )
    reactions = utility.read_file_table(
        path_file=path_file_reactions,
        names=[
            "identifier", "equation", "recon2m2", "metanetx",
            "enzyme_commission", "processes", "references"
        ],
        delimiter="\t"
    )
    # Compile information
    return {
        "genes": genes,
        "compartments": compartments,
        "metabolites": metabolites,
        "reactions": reactions
    }


def extract_compartments(compartments_source=None):
    """
    Extracts information from source about compartments

    arguments:
        compartments_source (list<dict>): source information about compartments

    returns:
        (dict<dict>): information about compartments

    raises:

    """

    compartments = {}
    for source_compartment in compartments_source:
        record = {
            "identifier": source_compartment["identifier"],
            "name": source_compartment["name"]
        }
        compartments[source_compartment["identifier"]] = record
    return compartments


def extract_processes(reactions_source=None):
    """
    Extracts information from source about processes

    arguments:
        reactions_source (list<dict>): source information about reactions

    returns:
        (dict<dict>): information about processes

    raises:

    """

    processes = {}
    for reaction_source in reactions_source:
        reaction_processes_names = extract_reaction_processes_names(
            reaction_source=reaction_source
        )
        for name in reaction_processes_names:
            # Determine whether a record exists for the process
            novelty = determine_process_name_novelty(
                name=name, processes=processes
            )
            if novelty:
                # Create and include a record for the process
                index = len(processes.keys())
                identifier = "P" + str(index + 1)
                record = {
                    "identifier": identifier,
                    "name": name
                }
                processes[identifier] = record
    return processes


def extract_reaction_processes_names(reaction_source=None):
    """
    Extracts names of a reaction's metabolic processes

    arguments:
        reaction_source (dict): source information about a reaction

    returns:
        (list<str>): names of a reaction's process

    raises:

    """

    # Separate references
    processes_source = reaction_source["processes"]
    processes = extract_reference_information(
        key="model:", references_source=processes_source
    )
    return processes


def extract_reference_information(key=None, references_source=None):
    """
    Extracts reference information

    arguments:
        key (str): identifier of reference information for a specific type
        references_source (str): source information about references

    returns:
        (list<str>): reference information for a specific type

    raises:

    """

    # Separate information for references
    references = references_source.split(";")
    # Filter identifiers for the reference
    pairs = list(filter(lambda pair: key in pair, references))
    # Remove key from identifiers
    identifiers = list(map(lambda pair: pair.replace(key, ""), pairs))
    # Return identifiers
    return identifiers


def determine_process_name_novelty(name=None, processes=None):
    """
    Determines whether a process' name is novel

    arguments:
        name (str): name of a process
        processes (dict<dict>): information about processes

    returns:
        (bool): whether process' name is novel

    raises:

    """

    for record in processes.values():
        if name == record["name"]:
            return False
    return True


def extract_metabolites(metabolites_source=None):
    """
    Extracts information from source about metabolites

    arguments:
        metabolites_source (list<dict>): source information about metabolites

    returns:
        (dict<dict>): information about metabolites

    raises:

    """

    metabolites = {}
    for metabolite_source in metabolites_source:
        record = extract_metabolite(metabolite_source=metabolite_source)
        metabolites[record["identifier"]] = record
    return metabolites


def extract_metabolite(metabolite_source=None):
    """
    Extracts information from source about a metabolite

    arguments:
        metabolite_source (dict): source information about a metabolite

    returns:
        (dict): information about a metabolite

    raises:

    """

    # Determine information
    identifier = metabolite_source["identifier"]
    name = metabolite_source["name"]
    formula = metabolite_source["formula"]
    mass = metabolite_source["mass"]
    charge = metabolite_source["charge"]
    references = extract_metabolite_references(
        identifier=identifier,
        references_source=metabolite_source["references"]
    )
    # Compile and return information
    return {
        "identifier": identifier,
        "name": name,
        "formula": formula,
        "mass": mass,
        "charge": charge,
        "references": references
    }


def extract_metabolite_references(identifier=None, references_source=None):
    """
    Extracts references from source about a metabolite

    arguments:
        identifier (str): identifier of metabolite in MetaNetX
        references_source (str): source information about a metabolite's
            references

    returns:
        (dict<str>): references about a metabolite

    raises:

    """

    # Collect identifiers for each reference
    chebi = extract_reference_information(
        references_source=references_source, key="chebi:"
    )
    bigg = extract_reference_information(
        references_source=references_source, key="bigg:"
    )
    metanetx_prior = extract_reference_information(
        references_source=references_source, key="deprecated:"
    )
    metanetx = [identifier] + metanetx_prior
    envipath = extract_reference_information(
        references_source=references_source, key="envipath:"
    )
    hmdb = extract_reference_information(
        references_source=references_source, key="hmdb:"
    )
    kegg = extract_reference_information(
        references_source=references_source, key="kegg:"
    )
    lipidmaps = extract_reference_information(
        references_source=references_source, key="lipidmaps:"
    )
    metacyc = extract_reference_information(
        references_source=references_source, key="metacyc:"
    )
    reactome = extract_reference_information(
        references_source=references_source, key="reactome:"
    )
    sabiork = extract_reference_information(
        references_source=references_source, key="sabiork:"
    )
    seed = extract_reference_information(
        references_source=references_source, key="seed:"
    )
    slm = extract_reference_information(
        references_source=references_source, key="slm:"
    )
    # Compile and return information
    return {
        "chebi": chebi,
        "bigg": bigg,
        "metanetx": metanetx,
        "envipath": envipath,
        "hmdb": hmdb,
        "kegg": kegg,
        "lipidmaps": lipidmaps,
        "metacyc": metacyc,
        "reactome": reactome,
        "sabiork": sabiork,
        "seed": seed,
        "slm": slm
    }


def extract_reactions(
    reactions_source=None, genes_source=None, processes=None
):
    """
    Extracts information from source about reactions

    arguments:
        reactions_source (list<dict>): source information about reactions
        genes_source (list<dict>): source information about genes
        processes (dict<dict>): information about processes

    returns:
        (dict<dict>): information about reactions

    raises:

    """

    reactions = {}
    for reaction_source in reactions_source:
        record = extract_reaction(
            reaction_source=reaction_source, genes_source=genes_source,
            processes=processes
        )
        reactions[record["identifier"]] = record
    return reactions


def extract_reaction(reaction_source=None, genes_source=None, processes=None):
    """
    Extracts information from source about a reaction

    arguments:
        reaction_source (dict): source information about a reaction
        genes_source (list<dict>): source information about genes
        processes (dict<dict>): information about processes

    returns:
        (dict): information about a reaction

    raises:

    """

    # Determine information
    identifier = reaction_source["identifier"]
    equation = reaction_source["equation"]
    reversibility = extract_reaction_reversibility(equation=equation)
    participants = extract_reaction_participants(equation=equation)
    processes = extract_reaction_processes(
        reaction_source=reaction_source,
        processes=processes
    )
    references = extract_reaction_references(
        recon2m2=reaction_source["recon2m2"],
        metanetx=reaction_source["metanetx"],
        enzyme_commission=reaction_source["enzyme_commission"],
        references_source=reaction_source["references"]
    )
    genes = extract_reaction_genes(
        identifier=identifier, genes_source=genes_source
    )
    # Compile and return information
    return {
        "identifier": identifier,
        "equation": equation,
        "reversibility": reversibility,
        "participants": participants,
        "processes": processes,
        "genes": genes,
        "references": references
    }


def extract_reaction_reversibility(equation=None):
    """
    Extracts information about a reaction's reversibility

    arguments:
        equation (str): a reaction's equation from MetaNetX

    returns:
        (bool): whether reaction is reversible

    raises:

    """

    if "<==>" in equation:
        return True
    else:
        return False


def extract_reaction_participants(equation=None):
    """
    Extracts information about a reaction's participants

    arguments:
        equation (str): a reaction's equation from MetaNetX

    returns:
        (list<dict>): information about a reaction's participants

    raises:

    """

    # Extract raw information about reaction's participants
    participants_raw = extract_reaction_equation_raw_participants_by_role(
        equation=equation
    )
    # Extract information about participants' role, coefficient, metabolite,
    # and compartment
    reactants = extract_reaction_participants_by_role(
        participants_raw=participants_raw["reactants"],
        role="reactant"
    )
    products = extract_reaction_participants_by_role(
        participants_raw=participants_raw["products"],
        role="product"
    )
    # Compile and return information
    return reactants + products


def extract_reaction_equation_raw_participants_by_role(equation=None):
    """
    Extracts raw information about a reaction's participants from its equation

    arguments:
        equation (str): a reaction's equation from MetaNetX

    returns:
        (dict<list<str>>): raw information about a reaction's participants by
            role

    raises:

    """

    # Separate information about participants' reactants from products
    # Determine reaction's directionality
    if "<==>" in equation:
        equation_sides = equation.split(" <==> ")
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]
    elif "-->" in equation:
        equation_sides = equation.split(" --> ")
        reactants_side = equation_sides[0]
        products_side = equation_sides[1]
    elif "<--" in equation:
        equation_sides = equation.split(" <-- ")
        reactants_side = equation_sides[1]
        products_side = equation_sides[0]
    # Separate information about individual participants
    reactants = reactants_side.split(" + ")
    products = products_side.split(" + ")
    # Compile and return information
    return {
        "reactants": reactants,
        "products": products
    }


def extract_reaction_participants_by_role(participants_raw=None, role=None):
    """
    Extracts information about a reaction's participants by role

    arguments:
        participants_raw (list<str>): raw information about a reaction's
            participants
        role (str): participants' role in reaction

    returns:
        (list<dict>): information about a reaction's participants

    raises:

    """

    # Extract information about participants
    participants = []
    for participant_raw in participants_raw:
        record = extract_reaction_participant_by_role(
            participant_raw=participant_raw, role=role
        )
        participants.append(record)
    return participants


def extract_reaction_participant_by_role(participant_raw=None, role=None):
    """
    Extracts information about a reaction's participant by role

    arguments:
        participant_raw (str): raw information about a reaction's participant
        role (str): participant's role in reaction

    returns:
        (dict): information about a reaction's participant

    raises:

    """

    # Separate information
    participant_split = participant_raw.split(" ")
    coefficient = float(participant_split[0])
    metabolite_compartment = participant_split[1].split("@")
    metabolite = metabolite_compartment[0]
    compartment = metabolite_compartment[1]
    # Compile and return information
    return {
        "metabolite": metabolite,
        "compartment": compartment,
        "coefficient": coefficient,
        "role": role
    }


def extract_reaction_genes(identifier=None, genes_source=None):
    """
    Extracts information about a reaction's genes

    arguments:
        identifier (str): identifier of a reaction
        genes_source (list<dict>): source information about genes

    returns:
        (list<str>): identifiers of a reaction's genes

    raises:

    """

    def match_reaction_gene(gene_record):
        return gene_record["reaction"] == identifier
    gene_source = utility.find(match_reaction_gene, genes_source)
    genes_references = gene_source["genes"]
    genes_split_one = genes_references.split(";")
    genes_split_two = genes_references.split("+")
    genes = []
    for gene_pair in genes_split_two:
        identifier = gene_pair.replace("gene:", "")
        genes.append(identifier)
    return genes


def extract_reaction_processes(reaction_source=None, processes=None):
    """
    Extracts identifiers of a reaction's metabolic processes

    arguments:
        reaction_source (dict): source information about a reaction
        processes (dict<dict>): information about processes

    returns:
        (list<str>): identifiers of a reaction's process

    raises:

    """

    reaction_processes_names = extract_reaction_processes_names(
        reaction_source=reaction_source
    )
    reaction_processes = []
    for name in reaction_processes_names:
        # Find process
        def match_reaction_process(record):
            return record["name"] == name
        process_record = utility.find(
            match_reaction_process, processes.values()
        )
        identifier = process_record["identifier"]
        reaction_processes.append(identifier)
    return reaction_processes


def extract_reaction_references(
    recon2m2=None, metanetx=None, enzyme_commission=None,
    references_source=None
):
    """
    Extracts references from source about a reaction

    arguments:
        recon2m2 (str): identifier of reaction in original version of model,
            Recon 2M.2
        metanetx (str): identifier of reaction in MetaNetX
        enzyme_commission (str): identifier of reaction in Enzyme Commission
        references_source (str): source information about a reaction's
            references

    returns:
        (dict<str>): references about a reaction

    raises:

    """

    # Collect identifiers for each reference
    rhea = extract_reference_information(
        references_source=references_source, key="rhea:"
    )
    bigg = extract_reference_information(
        references_source=references_source, key="bigg:"
    )
    metanetx_prior = extract_reference_information(
        references_source=references_source, key="deprecated:"
    )
    metanetx_current = [metanetx] + metanetx_prior
    kegg = extract_reference_information(
        references_source=references_source, key="kegg:"
    )
    metacyc = extract_reference_information(
        references_source=references_source, key="metacyc:"
    )
    reactome = extract_reference_information(
        references_source=references_source, key="reactome:"
    )
    sabiork = extract_reference_information(
        references_source=references_source, key="sabiork:"
    )
    seed = extract_reference_information(
        references_source=references_source, key="seed:"
    )
    # Compile and return information
    return {
        "recon2m2": recon2m2,
        "rhea": rhea,
        "bigg": bigg,
        "metanetx": metanetx_current,
        "enzyme_commission": enzyme_commission.split(";"),
        "kegg": kegg,
        "metacyc": metacyc,
        "reactome": reactome,
        "sabiork": sabiork,
        "seed": seed
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
    # Write information to file
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
    # Extract information about compartments
    compartments = extract_compartments(
        compartments_source=source["compartments"]
    )
    # Extract information about processes
    processes = extract_processes(reactions_source=source["reactions"])
    # Extract information about metabolites
    metabolites = extract_metabolites(metabolites_source=source["metabolites"])
    # Extract information about reactions
    reactions = extract_reactions(
        reactions_source=source["reactions"], genes_source=source["genes"],
        processes=processes
    )
    # Compile information
    metabolism_sets_entities = {
        "compartments": compartments,
        "processes": processes,
        "metabolites": metabolites,
        "reactions": reactions
    }
    #Write product information to file
    write_product(information=metabolism_sets_entities)


if __name__ == "__main__":
    main()
