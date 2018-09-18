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
            closeness[node] = 0.0
    for node in bottom:
        sp = dict(path_length(G, node))
        totsp = sum(sp.values())
        if totsp > 0.0 and len(G) > 1:
            closeness[node] = (n + 2 * (m - 1)) / totsp
            if normalized:
                s = (len(sp) - 1.0) / (len(G) - 1)
                closeness[node] *= s
        else:
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
    # Read information from file.
    with open(path_networkx, "rb") as file_source:
        network_elements = pickle.load(file_source)
    # Compile and return information.
    return {
        "network_elements": network_elements
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
    nodes_reactions=None,
    nodes_metabolites=None
):
    """
    Analyze individual nodes in a bipartite network.
    Analyze nodes for reactions and metabolites separately.

    arguments:
        network (object): instance of network in NetworkX
        nodes_reactions (list<str>): identifiers of reactions' nodes
        nodes_metabolites (list<str>): identifiers of metabolites' nodes

    raises:

    returns:
        (dict<dict<dict>>): information about network's nodes

    """

    # Reactions.
    reactions = analyze_network_nodes_group(
        network=network,
        nodes=nodes_reactions
    )
    # Metabolites.
    metabolites = analyze_network_nodes_group(
        network=network,
        nodes=nodes_metabolites
    )
    # Compile information.
    return {
        "reactions": reactions,
        "metabolites": metabolites
    }


def analyze_network_nodes_group(
    network=None,
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
    clustering coefficient

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about network's nodes

    """

    # Degree.
    degrees = determine_network_nodes_degrees(network=network, nodes=nodes)
    # Centrality.
    centralities = determine_network_nodes_centralities(
        network=network,
        nodes=nodes
    )
    # Cluster.
    clusters = determine_network_nodes_clusters(
        network=network,
        nodes=nodes
    )
    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "degree_in": degrees[node]["degree_in"],
            "degree_out": degrees[node]["degree_out"],
            "degree": degrees[node]["degree"],
            "centrality_degree": centralities[node]["degree"],
            "centrality_closeness": centralities[node]["closeness"],
            "centrality_betweenness": centralities[node]["betweenness"],
            "cluster_coefficient": clusters[node]["coefficient"]
        }
        collection[node] = entry
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
    centralities_closeness = (
        bipartite_closeness_centrality(
            network, nodes, normalized=True
        )
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
        nodes=nodes_reactions
    )
    # Metabolites.
    metabolites = analyze_bipartite_network_group(
        network=network,
        nodes=nodes_metabolites
    )
    # Compile information.
    return {
        "reactions": reactions,
        "metabolites": metabolites
    }


def analyze_bipartite_network_group(
    network=None,
    nodes=None
):
    """
    Analyze a bipartite group in a network.

    arguments:
        network (object): instance of network in NetworkX
        nodes (list<str>): identifiers of nodes in a bipartite set

    raises:

    returns:
        (dict<dict>): information about network

    """

    # Scale.

    #ntx.algorithms.bipartite.basic.density(network, nodes)
    #print(len(source["network_elements"]["nodes"]))
    #print(len(source["network_elements"]["nodes_reactions"]))
    #print(len(source["network_elements"]["nodes_metabolites"]))
    #print(network.number_of_nodes()) #network.order() = count of nodes
    #print(network.number_of_edges()) #network.size() = count of links
    # graph diameter (longest shortest path between any node pairs)

    # Connectivity.
    #list(ntx.connected_components(network))
    #degree_assortativity_coefficient()
    #ntx.algorithms.components.is_connected(network)
    #ntx.algorithms.components.is_strongly_connected(network)
    #ntx.algorithms.components.is_weakly_connected(network)
    #ntx.algorithms.components.connected_components(network)
    #ntx.algorithms.components.number_connected_components(network)
    # TODO: average node connectivity might not be accurate for bipartite
    #ntx.algorithms.connectivity.connectivity.average_node_connectivity(network)
    # TODO: I think all pairs node connectivity will work if I specify a single bipartite group
    # TODO: calculate average from these calculations... average for each bipartite group
    #ntx.algorithms.connectivity.connectivity.all_pairs_node_connectivity(network, nodes)

    # Centralization.
    # centralization (https://stackoverflow.com/questions/35243795/calculate-indegree-centralization-of-graph-with-python-networkx)
    # https://www.rdocumentation.org/packages/igraph/versions/0.7.1/topics/centralization
    # https://www.stat.washington.edu/~pdhoff/courses/567/Notes/l6_centrality.pdf

    # Cluster.
    # average clustering coefficient
    #ntx.clustering(network)
    #bipartite.average_clustering(network)
    #ntx.clustering(network, network.nodes())
    #ntx.average_clustering(network, network.nodes())
    if False:
        clusters = determine_network_nodes_clusters(
            network=network,
            nodes=nodes
        )
        network_indirect = network.to_undirected()
        # Need to use bipartite cluster algorithm.
        print("cluster bipartite for metabolites")
        cluster_metabolites = ntx.algorithms.bipartite.cluster.average_clustering(
            network, source["network_elements"]["nodes_metabolites"]
        )
        print(cluster_metabolites)
        print("cluster for metabolites")
        cluster = ntx.algorithms.cluster.average_clustering(
            network_indirect#, source["network_elements"]["nodes_metabolites"]
        )
        print(cluster)



    # Compile information.
    return {
        "reactions": reactions,
        "metabolites": metabolites
    }


    #average_degree_connectivity()

    # assortativity
    #print(ntx.node_connectivity(network))
    network_indirect = network.to_undirected()
    print(ntx.algorithms.bipartite.is_bipartite(network_indirect))
    print(ntx.algorithms.bipartite.is_bipartite_node_set(
        network_indirect, source["network_elements"]["nodes_reactions"]
    ))
    #print(ntx.all_pairs_node_connectivity(network))
    #print(ntx.is_connected(network)) # not implemented for directed graph
    #print(ntx.number_connected_components(network)) # not implemented for directed graph


    # Compile information.
    collection = {}
    for node in nodes:
        entry = {
            "identifier": node,
            "degree_in": degrees[node]["degree_in"],
            "degree_out": degrees[node]["degree_out"],
            "degree": degrees[node]["degree"],
            "centrality_degree": centralities[node]["degree"],
            "centrality_closeness": centralities[node]["closeness"],
            "centrality_betweenness": centralities[node]["betweenness"],
            "cluster_coefficient": clusters[node]["coefficient"]
        }
        collection[node] = entry
    return collection




# TODO: implement this function for centralization...

def calculate_network_centralization(
    type=None,
    bipartite=None,
    count_nodes=None,
    count_metabolites=None,
    count_reactions=None,
    centralities=None
):
    """
    Calculates the centralization of the entire network.

    Centralities need to be raw centralities, not standard centralities that
    are normalized to maxima.

    Normalization to the maximal centrality for the star graph depends on
    whether the graph is bipartite.
    The relevant star graph in a bipartite scenario would have n-1 nodes of the
    type opposite that of the central node.
    Same formulas will apply as a non-bipartite star graph, but must adjust the
    value of n to match the correct type.

    arguments:
        type (str): type of centrality, degree, closeness, betweenness,
            eigenvector
        count (int): count of nodes in network
        centralities (list<tuple<str,str,dict>>): centralities of all nodes in
            network

    raises:

    returns:
        (float): centralization of network

    """

    # TODO: re-calculate centralities of all nodes without normalization
    # TODO: non-normalized degree centrality is just the node degree...
    # TODO: non-normalized betweenness centrality is just
    # TODO: ntx.betweenness_centrality(G, normalized=False)

    pass


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
        nodes=source["network_elements"]["nodes"],
        links=source["network_elements"]["links"]
    )
    # Network is bipartite.
    # Store references to separate groups of nodes for reactions and
    # metabolites.
    # Analyze network's individual nodes.
    report_network_nodes = analyze_network_nodes(
        network=network,
        nodes_reactions=source["network_elements"]["nodes_reactions"],
        nodes_metabolites=source["network_elements"]["nodes_metabolites"]
    )
    print(list(report_network_nodes["metabolites"].values())[10])
    # Analyze entire network.
    if False:
        report_network = analyze_network(
            network=network,
            nodes_reactions=source["network_elements"]["nodes_reactions"],
            nodes_metabolites=source["network_elements"]["nodes_metabolites"]
        )

    # Network analysis.



    # Prepare reports.
    # TODO: Sort records within reports...


    if False:
        # Compile information.
        information = {
            "nodes_reactions": nodes_reactions,
            "nodes_metabolites": nodes_metabolites,
            "links": links
        }
        #Write product information to file.
        write_product(directory=directory, information=information)
