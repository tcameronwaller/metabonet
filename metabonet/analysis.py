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

#dir()
#importlib.reload()

###############################################################################
# Functionality


# TODO: this is a temporary bug fix for...
# TODO: ntx.algorithms.bipartite.centrality.closeness_centrality()

def bipartite_closeness_centrality(G, nodes, normalized=True):
    r"""Compute the closeness centrality for nodes in a bipartite network.

    The closeness of a node is the distance to all other nodes in the
    graph or in the case that the graph is not connected to all other nodes
    in the connected component containing that node.

    Parameters
    ----------
    G : graph
        A bipartite network

    nodes : list or container
        Container with all nodes in one bipartite node set.

    normalized : bool, optional
      If True (default) normalize by connected component size.

    Returns
    -------
    closeness : dictionary
        Dictionary keyed by node with bipartite closeness centrality
        as the value.

    See Also
    --------
    betweenness_centrality,
    degree_centrality
    sets,
    is_bipartite

    Notes
    -----
    The nodes input parameter must conatin all nodes in one bipartite node set,
    but the dictionary returned contains all nodes from both node sets.
    See :mod:`bipartite documentation <networkx.algorithms.bipartite>`
    for further details on how bipartite graphs are handled in NetworkX.


    Closeness centrality is normalized by the minimum distance possible.
    In the bipartite case the minimum distance for a node in one bipartite
    node set is 1 from all nodes in the other node set and 2 from all
    other nodes in its own set [1]_. Thus the closeness centrality
    for node `v`  in the two bipartite sets `U` with
    `n` nodes and `V` with `m` nodes is

    .. math::

        c_{v} = \frac{m + 2(n - 1)}{d}, \mbox{for} v \in U,

        c_{v} = \frac{n + 2(m - 1)}{d}, \mbox{for} v \in V,

    where `d` is the sum of the distances from `v` to all
    other nodes.

    Higher values of closeness  indicate higher centrality.

    As in the unipartite case, setting normalized=True causes the
    values to normalized further to n-1 / size(G)-1 where n is the
    number of nodes in the connected part of graph containing the
    node.  If the graph is not completely connected, this algorithm
    computes the closeness centrality for each connected part
    separately.

    References
    ----------
    .. [1] Borgatti, S.P. and Halgin, D. In press. "Analyzing Affiliation
        Networks". In Carrington, P. and Scott, J. (eds) The Sage Handbook
        of Social Network Analysis. Sage Publications.
        http://www.steveborgatti.com/research/publications/bhaffiliations.pdf
    """
    closeness = {}
    path_length = ntx.single_source_shortest_path_length
    top = set(nodes)
    bottom = set(G) - top
    n = float(len(top))
    m = float(len(bottom))
    for node in top:
        sp = dict(path_length(G, node))
        totsp = sum(sp.values())
        if totsp > 0.0 and len(G) > 1:
            closeness[node] = (m + 2 * (n - 1)) / totsp
            if normalized:
                s = (len(sp) - 1.0) / (len(G) - 1)
                closeness[node] *= s
        else:
            # TODO: correction of error in original
            closeness[node] = 0.0
    for node in bottom:
        sp = dict(path_length(G, node))
        totsp = sum(sp.values())
        if totsp > 0.0 and len(G) > 1:
            closeness[node] = (n + 2 * (m - 1)) / totsp
            if normalized:
                s = (len(sp) - 1.0) / (len(G) - 1)
                closeness[node] = (closeness[node] / s)
                closeness[node] *= s
        else:
            # TODO: correction of error in original
            closeness[node] = 0.0
    return closeness


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
    path_networkx = os.path.join(
        path_conversion, "network_elements_networkx.pickle"
    )
    #path_network = os.path.join(directory, "network")
    #path_nodes_reactions = os.path.join(path_network, "nodes_reactions.pickle")
    #path_nodes_metabolites = os.path.join(
    #    path_network, "nodes_metabolites.pickle"
    #)
    # Read information from file.
    with open(path_networkx, "rb") as file_source:
        information = pickle.load(file_source)
    # Compile and return information.
    return {
        "nodes": information["nodes"],
        "links": information["links"],
        "nodes_reactions_identifiers": (
            information["nodes_reactions_identifiers"]
        ),
        "nodes_reactions": information["nodes_reactions"],
        "nodes_metabolites_identifiers": (
            information["nodes_metabolites_identifiers"]
        ),
        "nodes_metabolites": information["nodes_metabolites"],
    }


def instantiate_networkx(nodes=None, links=None):
    """
    Converts information about network's nodes and links to format for
    NetworkX.

    arguments:
        nodes (list<tuple<str,dict>>): information about network's nodes
        links (list<tuple<str,str,dict>>): information about network's links

    raises:

    returns:
        (object): instance of network in NetworkX

    """

    network = ntx.DiGraph()
    # Method add_nodes_from accepts a list of tuples of node identifier (str),
    # node attributes (dict).
    network.add_nodes_from(nodes)
    # Method add_edges_from accepts a list of tuples of source node identifier
    # (str), target node identifier (str), link attributes (dict).
    network.add_edges_from(links)
    # Return reference.
    return network


def analyze_bipartite_network_nodes(
    network=None,
    nodes_reactions_identifiers=None,
    nodes_reactions=None,
    nodes_metabolites_identifiers=None,
    nodes_metabolites=None
):
    """
    Analyze individual nodes in a bipartite network.
    Analyze nodes for reactions and metabolites separately.

    arguments:
        network (object): instance of network in NetworkX
        nodes_reactions_identifiers (list<str>): identifiers of reactions'
            nodes
        nodes_reactions (dict<dict>): information about network's nodes for
            reactions
        nodes_metabolites_identifiers (list<str>): identifiers of metabolites'
            nodes
        nodes_metabolites (dict<dict>): information about network's nodes for
            metabolites

    raises:

    returns:
        (dict<dict<dict>>): information about network's nodes

    """

    # Reactions.
    reactions = analyze_network_nodes_group(
        network=network,
        nodes_identifiers=nodes_reactions_identifiers,
        nodes=nodes_reactions
    )
    # Metabolites.
    metabolites = analyze_network_nodes_group(
        network=network,
        nodes_identifiers=nodes_metabolites_identifiers,
        nodes=nodes_metabolites
    )
    # Compile information.
    return {
        "reactions": reactions,
        "metabolites": metabolites
    }


def analyze_network_nodes_group(
    network=None,
    nodes_identifiers=None,
    nodes=None
):
    """
    Analyze individual nodes from a bipartite group in a network.

    in degree
    out degree
    degree
    degree centrality
    closeness centrality
    betweenness centrality
    eigenvector centrality
    eccentricity centrality
    eccentricity
    clustering coefficient

    arguments:
        network (object): instance of network in NetworkX
        nodes_identifiers (list<str>): identifiers of nodes in a bipartite set
        nodes (dict<dict>): information about network's nodes in a bipartite
            set

    raises:

    returns:
        (dict<dict>): information about network's nodes

    """

    # Degree.
    degrees = determine_network_nodes_degrees(
        network=network, nodes=nodes_identifiers
    )
    # Centrality.
    centralities = determine_network_nodes_centralities(
        network=network,
        nodes=nodes_identifiers
    )
    # Distance.
    # Measurements of distance including diameter and radius are undefined for
    # non-connected or even less-connected networks.
    #distances = determine_network_nodes_distances(network=network, nodes=nodes)
    # Cluster.
    clusters = determine_network_nodes_clusters(
        network=network,
        nodes=nodes_identifiers
    )
    # Rank.
    # Compile information.
    collection = {}
    for node_identifier in nodes_identifiers:
        entry = {
            "identifier": node_identifier,
            "type": nodes[node_identifier]["type"],
            "entity": nodes[node_identifier]["entity"],
            "name": nodes[node_identifier]["name"],
            "degree_in": degrees[node_identifier]["degree_in"],
            "degree_out": degrees[node_identifier]["degree_out"],
            "degree": degrees[node_identifier]["degree"],
            "centrality_degree": centralities[node_identifier]["degree"],
            "centrality_closeness": centralities[node_identifier]["closeness"],
            "centrality_betweenness": centralities[node_identifier]["betweenness"],
            #"eccentricity": distances[node]["eccentricity"],
            "cluster_coefficient": clusters[node_identifier]["coefficient"]
        }
        collection[node_identifier] = entry
    return collection


def determine_network_nodes_degrees(
    network=None,
    nodes=None
):
    """
    Determine degrees of nodes in a network.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about degrees of network's nodes

    """

    # Determine degrees.
    degrees_in = dict(network.in_degree(nodes))
    degrees_out = dict(network.out_degree(nodes))
    degrees = dict(network.degree(nodes))
    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "degree_in": degrees_in[node],
            "degree_out": degrees_out[node],
            "degree": degrees[node]
        }
        collection[node] = entry
    return collection


def determine_network_nodes_centralities(
    network=None,
    nodes=None
):
    """
    Determine centralities of nodes in a network.

    degree centrality
    closeness centrality
    betweenness centrality
    eigenvector centrality (ambiguous definition for bipartite network)
    eccentricity centrality (implementation unavailable)

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about centralities of network's nodes

    """

    # Determine centralities.
    centralities_degree = (
        ntx.algorithms.bipartite.centrality.degree_centrality(network, nodes)
    )
    # Normalize centralities to the order of each node's component if network
    # is discontinuous.
    centralities_closeness = bipartite_closeness_centrality(
        network, nodes, normalized=True
    )
    centralities_betweenness = (
        ntx.algorithms.bipartite.centrality.betweenness_centrality(
            network, nodes
        )
    )
    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "degree": centralities_degree[node],
            "closeness": centralities_closeness[node],
            "betweenness": centralities_betweenness[node]
        }
        collection[node] = entry
    return collection


def determine_network_nodes_distances(
    network=None,
    nodes=None
):
    """
    Determine distances of nodes in a network.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about distances of network's nodes

    """

    # Determine distances.
    eccentricities = ntx.algorithms.distance_measures.eccentricity(
        network, nodes
    )
    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "eccentricity": eccentricities[node]
        }
        collection[node] = entry
    return collection


def determine_network_nodes_clusters(
    network=None,
    nodes=None
):
    """
    Determine cluster coefficients of nodes in a network.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about cluster coefficients of network's nodes

    """

    # Determine cluster coefficients.
    clusters = (
        ntx.algorithms.bipartite.cluster.clustering(
            network, nodes, mode="dot"
        )
    )
    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "coefficient": clusters[node]
        }
        collection[node] = entry
    return collection


def analyze_bipartite_network(
    network=None,
    nodes_reactions=None,
    nodes_metabolites=None
):
    """
    Analyze a bipartite network.

    arguments:
        network (object): instance of network in NetworkX
        nodes_reactions (list<str>): identifiers of reactions' nodes
        nodes_metabolites (list<str>): identifiers of metabolites' nodes

    raises:

    returns:
        (dict<dict>): information about network

    """

    # Reactions.
    reactions = analyze_bipartite_network_group(
        network=network,
        nodes_cis=nodes_reactions,
        nodes_trans=nodes_metabolites
    )
    # Metabolites.
    metabolites = analyze_bipartite_network_group(
        network=network,
        nodes_cis=nodes_metabolites,
        nodes_trans=nodes_reactions
    )
    # Compile information.
    return {
        "reactions": reactions,
        "metabolites": metabolites
    }


def analyze_bipartite_network_group(
    network=None,
    nodes_cis=None,
    nodes_trans=None
):
    """
    Analyze a bipartite network with primary relevance to a single bipartite
    set of nodes.

    arguments:
        network (object): instance of network in NetworkX
        nodes_cis (list<str>): identifiers of nodes in the bipartite set of
            current primary relevance
        nodes_trans (list<str>): identifiers of nodes in the bipartite set of
            current secondary relevance

    raises:

    returns:
        (dict<dict>): information about network

    """

    # Scale.
    scale = determine_network_scale(network=network, nodes=nodes_cis)
    # Connectivity.
    connectivity = determine_network_connectivity(
        network=network, nodes=nodes_cis
    )
    # Distance.
    # Measurements of distance including diameter and radius are undefined for
    # non-connected or even less-connected networks.
    #distance = determine_network_distance(network=network, nodes=nodes_cis)
    # Centralization.
    centralization = determine_bipartite_network_centralization(
        network=network,
        nodes_cis=nodes_cis,
        nodes_trans=nodes_trans,
        method="separate"
    )
    # Cluster.
    cluster_coefficient = ntx.algorithms.bipartite.cluster.average_clustering(
        network, nodes=nodes_cis, mode="dot"
    )
    # Modularity.
    # NetworkX might lack this functionality.
    if False:
        modules = (
            ntx.algorithms.community.modularity_max.greedy_modularity_communities(
                network
        ))
        modules_orders = []
        for module in sorted(modules, key=len, reverse=True):
            modules_orders.append(len(module))
    # Compile information.
    collection = {
        "order_total": scale["order_total"],
        "order_set": scale["order_set"],
        "size": scale["size"],
        "density": scale["density"],
        "connection": connectivity["connection"],
        "connection_strong": connectivity["connection_strong"],
        "connection_weak": connectivity["connection_weak"],
        "components_count": connectivity["components_count"],
        "components_orders": connectivity["components_orders"],
        "bipartite": connectivity["bipartite"],
        #"diameter": distance["diameter"],
        #"radius": distance["radius"],
        "centralization_degree": centralization["degree"],
        "centralization_closeness": centralization["closeness"],
        "centralization_betweenness": centralization["betweenness"],
        "cluster_coefficient": cluster_coefficient,
        #"modules_orders": modules_orders
    }
    return collection


def determine_network_scale(
    network=None,
    nodes=None
):
    """
    Determine metrics of a network's scale.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict): information about a network

    """

    # Determine scale.
    #network.number_of_nodes()
    #network.number_of_edges()
    order_total = network.order()
    order_set = len(nodes)
    size = network.size()
    density = ntx.algorithms.bipartite.basic.density(network, nodes)
    # Compile information.
    collection = {
        "order_total": order_total,
        "order_set": order_set,
        "size": size,
        "density": density
    }
    return collection


def determine_network_connectivity(
    network=None,
    nodes=None
):
    """
    Determine metrics of a network's connectivity.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict): information about a network

    """

    # Determine connectivity.
    # Assortativity requires specific algorithm for bipartite network.
    #ntx.algorithms.assortativity.degree_assortativity_coefficient()
    #ntx.algorithms.assortativity.average_degree_connectivity()
    # Average node connectivity requires specific algorithm for bipartite
    # network.
    #ntx.algorithms.connectivity.connectivity.average_node_connectivity(network)
    # Calculate all pairs node connectivity for each bipartite group, and
    # then determine average.
    # This algorithm is computationally expensive.
    #connectivities = (
    #    ntx.algorithms.connectivity.connectivity.all_pairs_node_connectivity(
    #        network, nodes
    #    )
    #)
    connection = (
        ntx.algorithms.components.is_connected(network.to_undirected())
    )
    connection_strong = (
        ntx.algorithms.components.is_strongly_connected(network)
    )
    connection_weak = ntx.algorithms.components.is_weakly_connected(network)
    components_count = (
        ntx.algorithms.components.number_connected_components(
            network.to_undirected()
        )
    )
    components = ntx.algorithms.components.connected_components(
        network.to_undirected()
    )
    components_orders = []
    for component in sorted(components, key=len, reverse=True):
        components_orders.append(len(component))
    bipartite = ntx.algorithms.bipartite.is_bipartite(network.to_undirected())
    # Compile information.
    collection = {
        "connection": connection,
        "connection_strong": connection_strong,
        "connection_weak": connection_weak,
        "components_count": components_count,
        "components_orders": components_orders,
        "bipartite": bipartite,
    }
    return collection


def determine_network_distance(
    network=None,
    nodes=None
):
    """
    Determine metrics of a network's distance.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict): information about a network

    """

    # Determine distance.
    #networkx.algorithms.distance_measures.radius

    diameter = ntx.algorithms.distance_measures.diameter(network)
    radius = ntx.algorithms.distance_measures.radius(network)
    # Compile information.
    collection = {
        "diameter": diameter,
        "radius": radius
    }
    return collection


def determine_bipartite_network_centralization(
    network=None,
    nodes_cis=None,
    nodes_trans=None,
    method=None
):
    """
    Calculates the centralization of a bipartite network.

    It is reasonable to calculate centralization of a bipartite network in
    multiple ways.
    1. Consider both bipartite sets of nodes together to create a single
    measurement of the network's centralization.
    2. Consider each bipartite set of nodes separately to create two
    measurements of the network's centralization relative to each set of nodes.

    arguments:
        network (object): instance of network in NetworkX
        nodes_cis (list<str>): identifiers of nodes in the bipartite set of
            current primary relevance
        nodes_trans (list<str>): identifiers of nodes in the bipartite set of
            current secondary relevance
        method (str): method to calcuate centralization of bipartite network,
            relevant to both bipartite sets of nodes together or separate

    raises:

    returns:
        (dict): information about a network

    """

    # Determine non-normal centralities.
    # Non-normal bipartite degree centrality is node degree.
    #degree = dict(network.degree(nodes_cis))
    # Non-normal bipartite closeness centrality is undefined.
    # Non-normal bipartite betweenness centrality omits final normalization.
    #betweenness = ntx.algorithms.centrality.betweenness_centrality(
    #    network, normalized=False
    #)

    # Determine normal centralities.
    # Formulas for centralization will be relevant to these normal
    # centralities.
    # Degree centrality.
    degree = ntx.algorithms.bipartite.centrality.degree_centrality(
        network, nodes_cis
    )
    # Closeness centrality.
    # Centralization formulas might not account for the supplemental
    # normalization to order of node's component.
    closeness = bipartite_closeness_centrality(
        network, nodes_cis, normalized=True
    )
    # Betweenness centrality.
    betweenness = ntx.algorithms.bipartite.centrality.betweenness_centrality(
        network, nodes_cis
    )
    # Calculate centralization.
    centralization = calculate_bipartite_network_centralization(
        centralities_degree=degree,
        centralities_closeness=closeness,
        centralities_betweenness=betweenness,
        nodes_cis=nodes_cis,
        nodes_trans=nodes_trans,
        method=method
    )
    # Compile information.
    collection = {
        "degree": centralization["degree"],
        "closeness": centralization["closeness"],
        "betweenness": centralization["betweenness"]
    }
    return collection


def calculate_bipartite_network_centralization(
    centralities_degree=None,
    centralities_closeness=None,
    centralities_betweenness=None,
    nodes_cis=None,
    nodes_trans=None,
    method=None
):
    """
    Calculates the centralization of a bipartite network relevant to both
    bipartite sets of nodes either together or separate.

    arguments:
        centralities_degree (dict<float>): degree centralities of nodes in
            network
        centralities_closeness (dict<float>): closeness centralities of nodes
            in network
        centralities_betweenness (dict<float>): betweenness centralities of
            nodes in network
        nodes_cis (list<str>): identifiers of nodes in the bipartite set of
            current primary relevance
        nodes_trans (list<str>): identifiers of nodes in the bipartite set of
            current secondary relevance
        method (str): method to calcuate centralization of bipartite network,
            relevant to both bipartite sets of nodes together or separate

    raises:

    returns:
        (dict<float>): measurements of a network's centralization

    """

    # Determine whether to measure centrality relative to bipartite sets of
    # nodes together or separate.
    if method == "together":
        centralization_degree = (
            calculate_bipartite_network_centralization_together(
                centralities=centralities_degree,
                method="degree",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
        centralization_closeness = (
            calculate_bipartite_network_centralization_together(
                centralities=centralities_closeness,
                method="closeness",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
        centralization_betweenness = (
            calculate_bipartite_network_centralization_together(
                centralities=centralities_betweenness,
                method="betweenness",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
    elif method == "separate":
        centralization_degree = (
            calculate_bipartite_network_centralization_separate(
                centralities=centralities_degree,
                method="degree",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
        centralization_closeness = (
            calculate_bipartite_network_centralization_separate(
                centralities=centralities_closeness,
                method="closeness",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
        centralization_betweenness = (
            calculate_bipartite_network_centralization_separate(
                centralities=centralities_betweenness,
                method="betweenness",
                nodes_cis=nodes_cis,
                nodes_trans=nodes_trans
            )
        )
    # Compile information.
    collection = {
        "degree": centralization_degree,
        "closeness": centralization_closeness,
        "betweenness": centralization_betweenness
    }
    return collection


def calculate_bipartite_network_centralization_together(
    centralities=None,
    method=None,
    nodes_cis=None,
    nodes_trans=None
):
    """
    Calculates the centralization of a bipartite network relevant to both
    bipartite sets of nodes together.

    This procedure calculates the double mode centralization for both bipartite
    sets of nodes within a bipartite network.
    Formulas use normal centralities for these nodes.

    arguments:
        centralities (dict<float>): centralities of nodes in network
        method (str): method to measure centralities, degree, closeness, or
            betweenness
        nodes_cis (list<str>): identifiers of nodes in the bipartite set of
            current primary relevance
        nodes_trans (list<str>): identifiers of nodes in the bipartite set of
            current secondary relevance

    raises:

    returns:
        (float): centralization of network

    """

    # Define formulas to calculate and normalize centralizations relevant to
    # both bipartite sets of nodes together.
    # Centralities have normalization for bipartite networks.
    def calculate_centralization(
        centralities=None, nodes_cis=None, nodes_trans=None
    ):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        maximum = max(list(centralities.values()))
        n = no + ni
        return (n * maximum) - sum(list(centralities.values()))
    def maximize_degree(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        return (
            (no+ni-1)-
            ((no-1)/ni)-
            ((ni+no-1)/no)
        )
    def maximize_closeness(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        if no > ni:
            return (
                (2*(ni-1)*((ni+no-2)/((3*ni)+(4*no)-8)))+
                (2*(no-ni)*(((2*ni)-1)/((5*ni)+(2*no)-6)))+
                (2*(ni-1)*((no-2)/((2*ni)+(3*no)-6)))+
                (2*((ni-1)/(no+(4*ni)-4)))
            )
        else:
            return (
                (2*(no-1)*((ni+no-4)/((3*ni)+(4*no)-8)))+
                (2*(no-1)*((no-2)/((2*ni)+(3*no)-6)))+
                (2*(no-1)*((ni-no+1)/((2*ni)+(3*no)-4)))
            )
    def maximize_betweenness(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        if no > ni:
            return (
                (no+ni-1)-
                (((ni-1)*(no+ni-2)+(0.5*(no-ni)*(no+(3*ni)-3)))/
                ((0.5*no*(no-1))+(0.5*(ni-1)*(ni-2))+((no-1)*(ni-1))))
            )
        else:
            return (
                (no+ni-1)-
                (((no-1)*(no+ni-2))/(2*(no-1)*ni-1))
            )
    # Calculate and normalize centralizations.
    numerator = calculate_centralization(
        centralities=centralities,
        nodes_cis=nodes_cis,
        nodes_trans=nodes_trans
    )
    if method == "degree":
        denominator = maximize_degree(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    elif method == "closeness":
        denominator = maximize_closeness(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    elif method == "betweenness":
        denominator = maximize_betweenness(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    centralization = numerator / denominator
    # Return information.
    return centralization


def calculate_bipartite_network_centralization_separate(
    centralities=None,
    method=None,
    nodes_cis=None,
    nodes_trans=None
):
    """
    Calculates the centralization of a bipartite network relevant to both
    bipartite sets of nodes separate.

    This procedure calculates the single mode centralization for a single
    bipartite set of nodes within a bipartite network.
    Formulas use normal centralities for these nodes.

    arguments:
        centralities (dict<float>): centralities of nodes in network
        method (str): method to measure centralities, degree, closeness, or
            betweenness
        nodes_cis (list<str>): identifiers of nodes in the bipartite set of
            current primary relevance
        nodes_trans (list<str>): identifiers of nodes in the bipartite set of
            current secondary relevance

    raises:

    returns:
        (float): centralization of network

    """

    # Define formulas to calculate and normalize centralizations relevant to
    # both bipartite sets of nodes separate.
    # Centralities have normalization for bipartite networks.
    def calculate_centralization(
        centralities=None, nodes_cis=None, nodes_trans=None
    ):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        # Filter centralities to consider only nodes in cis bipartite set.
        centralities_cis = {}
        for node in centralities.keys():
            if node in nodes_cis:
                centralities_cis[node] = centralities[node]
        maximum = max(list(centralities_cis.values()))
        centralization = (no * maximum) - sum(list(centralities_cis.values()))
        return centralization
    def maximize_degree(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        return ((ni-1)*(no-1))
    def maximize_closeness(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        if no > ni:
            return (
                (((ni-1)*(no-2))/((2*no)-3))+
                (((ni-1)*(no-ni))/(no+ni-2))
            )
        else:
            return (
                ((no-2)*(no-1))/((2*no)-3)
            )
    def maximize_betweenness(nodes_cis=None, nodes_trans=None):
        no = len(nodes_cis)
        ni = len(nodes_trans)
        if no > ni:
            return (
                (2*((no-1)^2)*(ni-1))
            )
        else:
            return (
                (no-1)*
                (
                    (0.5*ni*(ni-1))+
                    (0.5*(no-1)*(no-2))+
                    ((no-1)*(ni-1))
                )
            )
    # Calculate and normalize centralizations.
    numerator = calculate_centralization(
        centralities=centralities,
        nodes_cis=nodes_cis,
        nodes_trans=nodes_trans
    )
    if method == "degree":
        denominator = maximize_degree(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    elif method == "closeness":
        denominator = maximize_closeness(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    elif method == "betweenness":
        denominator = maximize_betweenness(
            nodes_cis=nodes_cis,
            nodes_trans=nodes_trans
        )
    centralization = numerator / denominator
    # Return information.
    return centralization


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
    path = os.path.join(directory, "analysis")
    utility.confirm_path_directory(path)
    path_nodes_reactions = os.path.join(path, "nodes_reactions.tsv")
    path_nodes_metabolites = os.path.join(path, "nodes_metabolites.tsv")
    path_network_reactions = os.path.join(path, "network_reactions.tsv")
    path_network_metabolites = os.path.join(path, "network_metabolites.tsv")
    # Write information to file.
    utility.write_file_table(
        information=information["report_nodes_reactions"],
        path_file=path_nodes_reactions,
        names=information["report_nodes_reactions"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["report_nodes_metabolites"],
        path_file=path_nodes_metabolites,
        names=information["report_nodes_metabolites"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["report_network_reactions"],
        path_file=path_network_reactions,
        names=information["report_network_reactions"][0].keys(),
        delimiter="\t"
    )
    utility.write_file_table(
        information=information["report_network_metabolites"],
        path_file=path_network_metabolites,
        names=information["report_network_metabolites"][0].keys(),
        delimiter="\t"
    )


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to define network's elements.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Instantiate network in NetworkX.
    network = instantiate_networkx(
        nodes=source["nodes"],
        links=source["links"]
    )
    # Network is bipartite.
    # Store references to separate groups of nodes for reactions and
    # metabolites.
    # Analyze network's individual nodes.
    report_nodes = analyze_bipartite_network_nodes(
        network=network,
        nodes_reactions_identifiers=source["nodes_reactions_identifiers"],
        nodes_reactions=source["nodes_reactions"],
        nodes_metabolites_identifiers=source["nodes_metabolites_identifiers"],
        nodes_metabolites=source["nodes_metabolites"]
    )
    #print(list(report_nodes["metabolites"].values())[10])
    # Analyze entire network.
    report_network = analyze_bipartite_network(
        network=network,
        nodes_reactions=source["nodes_reactions_identifiers"],
        nodes_metabolites=source["nodes_metabolites_identifiers"]
    )
    #print(report_network)
    # Prepare reports.
    # TODO: Sort records within reports...
    # Compile information.
    information = {
        "report_nodes_reactions": list(report_nodes["reactions"].values()),
        "report_nodes_metabolites": list(
            report_nodes["metabolites"].values()
        ),
        "report_network_reactions": [report_network["reactions"]],
        "report_network_metabolites": [report_network["metabolites"]]
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
