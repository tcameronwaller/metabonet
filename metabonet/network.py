"""
Evaluate candidacy of metabolites and reactions for representation in network.

Candidate entities, metabolites and reactions, are entities that are elligible
candidates for representation in the network.
An entity's candidacy depends on the context of interest.
This context includes relevance compartmentalization in general, relevance of
specific compartments and processes, and relevance of individual entities.
An entity's candidacy also depends on the candidacies of other entities to
which the entity relates.
A reaction's candidacy depends on the candidacies of metabolites that
participate in it.
A metabolite's candidacy depends on the candidacies of reactions in which
it participates.

Title:

    candidacy

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
    argparse: Package to interpret user parameters from terminal.
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
import os
import pickle
import copy

# Relevant
import networkx as ntx

# Custom
import utility
import conversion
import analysis

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
    path_conversion = os.path.join(directory, "conversion")
    path_compartments = os.path.join(path_conversion, "compartments.pickle")
    path_processes = os.path.join(path_conversion, "processes.pickle")
    path_reactions = os.path.join(path_conversion, "reactions.pickle")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    path_candidacy = os.path.join(directory, "candidacy")
    path_reactions_candidacy = os.path.join(path_candidacy, "reactions.pickle")
    path_metabolites_candidacy = os.path.join(
        path_candidacy, "metabolites.pickle"
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
    with open(path_reactions_candidacy, "rb") as file_source:
        reactions_candidacy = pickle.load(file_source)
    with open(path_metabolites_candidacy, "rb") as file_source:
        metabolites_candidacy = pickle.load(file_source)
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites,
        "reactions_candidacy": reactions_candidacy,
        "metabolites_candidacy": metabolites_candidacy,
    }


def collect_reactions_nodes(
    reactions_candidacy=None,
    reactions=None,
    metabolites=None,
    compartments=None,
    processes=None
):
    """
    Collects information about reactions' nodes.

    arguments:
        reactions_candidacy (dict<dict>): information about candidate reactions
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes

    raises:

    returns:
        (dict<dict>): information about reactions' nodes

    """

    reactions_nodes = {}
    for reaction_candidacy in reactions_candidacy.values():
        reaction_node = define_reaction_node(
            reaction_candidacy_identifier=reaction_candidacy["identifier"],
            reactions_candidacy=reactions_candidacy,
            reactions=reactions,
            metabolites=metabolites,
            compartments=compartments,
            processes=processes
        )
        reactions_nodes[reaction_node["identifier"]] = reaction_node
    return reactions_nodes


def define_reaction_node(
    reaction_candidacy_identifier=None,
    reactions_candidacy=None,
    reactions=None,
    metabolites=None,
    compartments=None,
    processes=None
):
    """
    Defines information about a reaction's node.

    arguments:
        reaction_candidacy_identifier (str): identifier of a candidate reaction
        reactions_candidacy (dict<dict>): information about candidate reactions
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes

    raises:

    returns:
        (dict): information about a reaction's node

    """

    # Access information.
    reaction_candidacy = reactions_candidacy[reaction_candidacy_identifier]
    reaction = reactions[reaction_candidacy["reaction"]]
    processes_names = utility.collect_values_from_records_in_reference(
        key="name",
        identifiers=reaction["processes"],
        reference=processes
    )
    # Compile information.
    reaction_node = {
        "identifier": reaction_candidacy["identifier"],
        "entity": reaction["identifier"],
        "name": reaction_candidacy["name"],
        "reversibility": reaction_candidacy["reversibility"],
        "replicates": reaction_candidacy["replicates"],
        "processes": processes_names,
        "type": "reaction"
    }
    # Return information.
    return reaction_node


def collect_metabolites_nodes(
    metabolites_candidacy=None,
    reactions=None,
    metabolites=None,
    compartments=None,
    processes=None
):
    """
    Collects information about metabolites' nodes.

    arguments:
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes

    raises:

    returns:
        (dict<dict>): information about metabolites' nodes

    """

    metabolites_nodes = {}
    for metabolite_candidacy in metabolites_candidacy.values():
        metabolite_node = define_metabolite_node(
            metabolite_candidacy_identifier=metabolite_candidacy["identifier"],
            metabolites_candidacy=metabolites_candidacy,
            reactions=reactions,
            metabolites=metabolites,
            compartments=compartments,
            processes=processes
        )
        metabolites_nodes[metabolite_node["identifier"]] = metabolite_node
    return metabolites_nodes


def define_metabolite_node(
    metabolite_candidacy_identifier=None,
    metabolites_candidacy=None,
    reactions=None,
    metabolites=None,
    compartments=None,
    processes=None
):
    """
    Defines information about a metabolite's node.

    arguments:
        metabolite_candidacy_identifier (str): identifier of a candidate
            metabolite
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes

    raises:

    returns:
        (dict): information about a metabolite's node

    """

    # Access information.
    metabolite_candidacy = (
        metabolites_candidacy[metabolite_candidacy_identifier]
    )
    metabolite = metabolites[metabolite_candidacy["metabolite"]]
    compartment = metabolite_candidacy["compartment"]
    if compartment != "null":
        compartment_name = compartments[compartment]
    else:
        compartment_name = "null"
    # Compile information.
    metabolite_node = {
        "identifier": metabolite_candidacy["identifier"],
        "entity": metabolite["identifier"],
        "name": metabolite_candidacy["name"],
        "compartment": compartment_name,
        "formula": metabolite["formula"],
        "mass": metabolite["mass"],
        "charge": metabolite["charge"],
        "reference_hmdb": metabolite["references"]["hmdb"],
        "reference_pubchem": metabolite["references"]["pubchem"],
        "replication": metabolite_candidacy["replication"],
        "type": "metabolite"
    }
    # Return information.
    return metabolite_node


def collect_links(
    reactions_candidacy=None,
    metabolites_candidacy=None
):
    """
    Collects information about reactions' nodes.

    arguments:
        reactions_candidacy (dict<dict>): information about candidate reactions
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites

    raises:

    returns:
        (dict<dict>): information about links between reactions and metabolites

    """

    links = {}
    for reaction_candidacy in reactions_candidacy.values():
        links_reaction = define_reaction_links(
            reaction_candidacy_identifier=reaction_candidacy["identifier"],
            reactions_candidacy=reactions_candidacy,
            metabolites_candidacy=metabolites_candidacy
        )
        # Include links in collection.
        for link in links_reaction:
            candidacy = determine_link_candidacy(
                link=link,
                links=links,
                reactions_candidacy=reactions_candidacy,
                metabolites_candidacy=metabolites_candidacy
            )
            if candidacy:
                links[link["identifier"]] = link
    return links


def define_reaction_links(
    reaction_candidacy_identifier=None,
    reactions_candidacy=None,
    metabolites_candidacy=None
):
    """
    Defines information about a reaction's links.

    arguments:
        reaction_candidacy_identifier (str): identifier of a candidate reaction
        reactions_candidacy (dict<dict>): information about candidate reactions
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites

    raises:

    returns:
        (list<dict>): information about reaction's links

    """

    # Access information.
    reaction_candidacy = reactions_candidacy[reaction_candidacy_identifier]
    participants = reaction_candidacy["participants"]
    # Create links for reaction's participants.
    links_reaction = []
    for participant in participants:
        links_participant = define_reaction_participant_links(
            participant=participant,
            reversibility=reaction_candidacy["reversibility"],
            reaction_candidacy_identifier=reaction_candidacy_identifier
        )
        links_reaction.extend(links_participant)
    return links_reaction


def define_reaction_participant_links(
    participant=None,
    reversibility=None,
    reaction_candidacy_identifier=None
):
    """
    Defines information about a reaction's participant's links.

    arguments:
        participant (dict<str>): information about a metabolite and a
            compartment that participate in a reaction
        reversibility (bool): whether a reaction is reversible
        reaction_candidacy_identifier (str): identifier of a candidate reaction

    raises:

    returns:
        (list<dict>): information about reaction's participant's links

    """

    # Access information.
    metabolite_candidacy_identifier = participant["metabolite_candidacy"]
    replication = participant["replication"]
    role = participant["role"]
    # Determine which links to define.
    if reversibility:
        # Reaction is reversible.
        # Define links for participant as both reactant and product.
        link_forward = define_link(
            source=metabolite_candidacy_identifier,
            target=reaction_candidacy_identifier,
            role=role,
            replication=replication
        )
        link_reverse = define_link(
            source=reaction_candidacy_identifier,
            target=metabolite_candidacy_identifier,
            role=role,
            replication=replication
        )
        links_participant = [link_forward, link_reverse]
    else:
        # Reaction is irreversible.
        # Define link for participant according to its role as reactant or
        # product.
        if role == "reactant":
            link = define_link(
                source=metabolite_candidacy_identifier,
                target=reaction_candidacy_identifier,
                role=role,
                replication=replication
            )
        elif role == "product":
            link = define_link(
                source=reaction_candidacy_identifier,
                target=metabolite_candidacy_identifier,
                role=role,
                replication=replication
            )
        links_participant = [link]
    return links_participant


def define_link(
    source=None,
    target=None,
    role=None,
    replication=None
):
    """
    Defines information about a link.

    arguments:
        source (str): identifier of a node that is link's source
        target (str): identifier of a node that is link's target
        role (str): role of link's participant, either reactant or product
        replication (bool): whether candidate metabolite has simplification by
            replication

    raises:

    returns:
        (dict<str>): information about a link

    """

    identifier = source + "_-_" + target
    link = {
        "identifier": identifier,
        "source": source,
        "target": target,
        "role": role,
        "replication": replication
    }
    return link


def determine_link_candidacy(
    link=None,
    links=None,
    reactions_candidacy=None,
    metabolites_candidacy=None
):
    """
    Determines whether link is a candidate.

    arguments:
        link (dict<str>): information about a link
        links (dict<dict>): information about links between reactions and
            metabolites
        reactions_candidacy (dict<dict>): information about candidate reactions
        metabolites_candidacy (dict<dict>): information about candidate
            metabolites

    raises:

    returns:
        (bool): whether link is a candidate

    """

    # Determine whether link refers to valid nodes.
    source_reactions = link["source"] in reactions_candidacy.keys()
    source_metabolites = link["source"] in metabolites_candidacy.keys()
    source = source_reactions or source_metabolites
    target_reactions = link["target"] in reactions_candidacy.keys()
    target_metabolites = link["target"] in metabolites_candidacy.keys()
    target = target_reactions or target_metabolites
    # Determine whether link is novel in collection.
    novelty = link["identifier"] not in links.keys()
    # Determine whether link is a candidate.
    candidacy = source and target and novelty
    return candidacy


def select_network_component_elements(
    component=None,
    nodes_reactions=None,
    nodes_metabolites=None,
    links=None,
):
    """
    Determines whether to select the main component of a network.

    arguments:
        component (bool): whether to select network's main component
        nodes_reactions (dict<dict>): information about reactions' nodes
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        links (dict<dict>): information about links between nodes for reactions
            and metabolites

    raises:

    returns:
        (dict<dict<dict>>): information about network's nodes and links

    """

    if component:
        # Convert information format for export to NetworkX.
        networkx = conversion.convert_networkx(
            nodes_reactions=nodes_reactions,
            nodes_metabolites=nodes_metabolites,
            links=links
        )
        # Instantiate network in NetworkX.
        network = analysis.instantiate_networkx(
            nodes=networkx["nodes"],
            links=networkx["links"]
        )
        # Select network's main component.
        components = ntx.algorithms.components.connected_components(
            network.to_undirected()
        )
        components_sort = sorted(components, key=len, reverse=True)
        component_main = components_sort[0]
        # Induce subnetwork for component.
        subnetwork = ntx.DiGraph.subgraph(network, list(component_main))
        # Extract identifiers of nodes and links from subnetwork.
        nodes_identifiers = []
        for node, data in subnetwork.nodes.items():
            nodes_identifiers.append(data["identifier"])
        links_identifiers = []
        for link, data in subnetwork.edges.items():
            links_identifiers.append(data["identifier"])
        # Filter nodes and links by identifiers.
        nodes_reactions_component = utility.filter_entries_identifiers(
            identifiers=nodes_identifiers,
            entries_original=nodes_reactions
        )
        nodes_metabolites_component = utility.filter_entries_identifiers(
            identifiers=nodes_identifiers,
            entries_original=nodes_metabolites
        )
        links_component = utility.filter_entries_identifiers(
            identifiers=links_identifiers,
            entries_original=links
        )
        # Compile and return information.
        return {
            "nodes_reactions": nodes_reactions_component,
            "nodes_metabolites": nodes_metabolites_component,
            "links": links_component
        }
    else:
        # Compile and return information.
        return {
            "nodes_reactions": nodes_reactions,
            "nodes_metabolites": nodes_metabolites,
            "links": links
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
    path = os.path.join(directory, "network")
    utility.confirm_path_directory(path)
    path_nodes_reactions = os.path.join(path, "nodes_reactions.pickle")
    path_nodes_metabolites = os.path.join(path, "nodes_metabolites.pickle")
    path_links = os.path.join(path, "links.pickle")
    # Write information to file.
    with open(path_nodes_reactions, "wb") as file_product:
        pickle.dump(information["nodes_reactions"], file_product)
    with open(path_nodes_metabolites, "wb") as file_product:
        pickle.dump(information["nodes_metabolites"], file_product)
    with open(path_links, "wb") as file_product:
        pickle.dump(information["links"], file_product)


###############################################################################
# Procedure


def execute_procedure(
    component=None,
    directory=None
):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to define network's elements.

    arguments:
        component (bool): whether to select network's main component
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Define network's nodes for reactions.
    nodes_reactions = collect_reactions_nodes(
        reactions_candidacy=source["reactions_candidacy"],
        reactions=source["reactions"],
        metabolites=source["metabolites"],
        compartments=source["compartments"],
        processes=source["processes"]
    )
    # Define network's nodes for metabolites.
    nodes_metabolites = collect_metabolites_nodes(
        metabolites_candidacy=source["metabolites_candidacy"],
        reactions=source["reactions"],
        metabolites=source["metabolites"],
        compartments=source["compartments"],
        processes=source["processes"]
    )
    # Define network's links.
    links = collect_links(
        reactions_candidacy=source["reactions_candidacy"],
        metabolites_candidacy=source["metabolites_candidacy"]
    )
    # Determine whether to select network's main component.
    component_elements = select_network_component_elements(
        component=component,
        nodes_reactions=nodes_reactions,
        nodes_metabolites=nodes_metabolites,
        links=links,
    )
    # Compile information.
    information = {
        "nodes_reactions": component_elements["nodes_reactions"],
        "nodes_metabolites": component_elements["nodes_metabolites"],
        "links": component_elements["links"]
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
