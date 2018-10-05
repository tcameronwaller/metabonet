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
import math

# Relevant.
import numpy
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pyplot
import wordcloud

# Custom
import utility

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
    path_analysis = os.path.join(directory, "analysis")
    path_measurement = os.path.join(directory, "measurement")
    # Nodes.
    path_nodes_compartments_yes_hubs_yes = os.path.join(
        path_analysis, "nodes_metabolites_compartments-yes_hubs-yes.pickle"
    )
    path_nodes_compartments_yes_hubs_no = os.path.join(
        path_analysis, "nodes_metabolites_compartments-yes_hubs-no.pickle"
    )
    path_nodes_compartments_no_hubs_yes = os.path.join(
        path_analysis, "nodes_metabolites_compartments-no_hubs-yes.pickle"
    )
    path_nodes_compartments_no_hubs_no = os.path.join(
        path_analysis, "nodes_metabolites_compartments-no_hubs-no.pickle"
    )
    # Measurements.
    path_measurements_one = os.path.join(path_measurement, "study_one.pickle")
    path_measurements_two = os.path.join(path_measurement, "study_two.pickle")
    path_measurements_three = os.path.join(
        path_measurement, "study_three.pickle"
    )
    path_measurements_four = os.path.join(
        path_measurement, "study_four.pickle"
    )
    path_measurements_five = os.path.join(
        path_measurement, "study_five.pickle"
    )
    # Read information from file.
    # Nodes.
    with open(path_nodes_compartments_yes_hubs_yes, "rb") as file_source:
        nodes_compartments_yes_hubs_yes = pickle.load(file_source)
    with open(path_nodes_compartments_yes_hubs_no, "rb") as file_source:
        nodes_compartments_yes_hubs_no = pickle.load(file_source)
    with open(path_nodes_compartments_no_hubs_yes, "rb") as file_source:
        nodes_compartments_no_hubs_yes = pickle.load(file_source)
    with open(path_nodes_compartments_no_hubs_no, "rb") as file_source:
        nodes_compartments_no_hubs_no = pickle.load(file_source)
    # Measurements.
    with open(path_measurements_one, "rb") as file_source:
        measurements_one = pickle.load(file_source)
    with open(path_measurements_two, "rb") as file_source:
        measurements_two = pickle.load(file_source)
    with open(path_measurements_three, "rb") as file_source:
        measurements_three = pickle.load(file_source)
    with open(path_measurements_four, "rb") as file_source:
        measurements_four = pickle.load(file_source)
    with open(path_measurements_five, "rb") as file_source:
        measurements_five = pickle.load(file_source)
    # Compile and return information.
    return {
        "nodes_compartments_yes_hubs_yes": nodes_compartments_yes_hubs_yes,
        "nodes_compartments_yes_hubs_no": nodes_compartments_yes_hubs_no,
        "nodes_compartments_no_hubs_yes": nodes_compartments_no_hubs_yes,
        "nodes_compartments_no_hubs_no": nodes_compartments_no_hubs_no,
        "measurements_one": measurements_one,
        "measurements_two": measurements_two,
        "measurements_three": measurements_three,
        "measurements_four": measurements_four,
        "measurements_five": measurements_five
    }


def define_font_properties():
    """
    Defines font properties.

    arguments:

    raises:

    returns:
        (dict<object>): references to definitions of font properties

    """

    font_one = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=1000,
        size=30
    )
    font_two = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=1000,
        size=25
    )
    font_three = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=1000,
        size=20
    )
    # Compile and return references.
    return {
        "font_one": font_one,
        "font_two": font_two,
        "font_three": font_three
    }


def define_color_properties():
    """
    Defines color properties.

    arguments:

    raises:

    returns:
        (dict<tuple>): references to definitions of color properties

    """

    # Black.
    black = (0.0, 0.0, 0.0, 1.0)
    # White.
    white = (1.0, 1.0, 1.0, 1.0)
    white_faint = (1.0, 1.0, 1.0, 0.75)
    # Blue.
    blue = (0.0, 0.2, 0.5, 1.0)
    blue_faint = (0.0, 0.2, 0.5, 0.75)
    # Orange.
    orange = (1.0, 0.6, 0.2, 1.0)
    orange_faint = (1.0, 0.6, 0.2, 0.75)
    # Compile and return references.
    return {
        "black": black,
        "white": white,
        "white_faint": white_faint,
        "blue": blue,
        "blue_faint": blue_faint,
        "orange": orange,
        "orange_faint": orange_faint
    }


def plot_degrees(
    nodes_hubs_yes=None,
    nodes_hubs_no=None,
    fonts=None,
    colors=None
):
    """
    Creates a chart to represent distributions of degrees of metabolites'
    nodes.

    arguments:
        nodes_hubs_yes (dict<dict>): information about metabolites' nodes
        nodes_hubs_no (dict<dict>): information about metabolites' nodes
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    series_one = utility.collect_value_from_records(
        key="degree",
        records=nodes_hubs_yes.values()
    )
    series_two = utility.collect_value_from_records(
        key="degree",
        records=nodes_hubs_no.values()
    )
    chart_degrees = plot_two_distributions_histograms(
        series_one=series_one,
        series_two=series_two,
        name_one="+ Hubs",
        name_two="-- Hubs",
        fonts=fonts,
        colors=colors
    )
    # Return reference to figure.
    return chart_degrees


def plot_two_distributions_histograms(
    series_one=None,
    series_two=None,
    name_one=None,
    name_two=None,
    fonts=None,
    colors=None
):
    """
    Creates a histogram chart to represent frequency distributions of two
    series.

    arguments:
        series_one (list<int>): series of counts
        series_two (list<int>): series of counts
        name_one (str): name of distribution
        name_two (str): name of distribution
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    # Determine histogram bins.
    hist, bin_edges = numpy.histogram(series_one, bins=50)
    # Create chart.
    figure = pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = pyplot.axes()
    values_one, bins_one, patches_one = axes.hist(
        series_one,
        bins=bin_edges,
        histtype="bar",
        align="left",
        orientation="vertical",
        rwidth=0.35,
        log=True,
        color=colors["blue"],
        label=name_one,
        stacked=False
    )
    values_two, bins_two, patches_two = axes.hist(
        series_two,
        bins=bin_edges,
        histtype="bar",
        align="mid",
        orientation="vertical",
        rwidth=0.35,
        log=True,
        color=colors["orange"],
        label=name_two,
        stacked=False
    )
    axes.legend(
        loc="upper right",
        markerscale=2.5,
        markerfirst=True,
        prop=fonts["font_one"],
        edgecolor=colors["black"]
    )
    axes.set_xlabel(
        xlabel="Node Degree",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["font_one"]
    )
    axes.set_ylabel(
        ylabel="Count of Nodes",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["font_one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["font_two"].get_size(),
        labelcolor=colors["black"]
    )
    return figure


def plot_ranks(
    color=None,
    count=None,
    nodes_metabolites=None,
    fonts=None,
    colors=None
):
    """
    Creates a chart to represent ranks of metabolites' nodes.

    arguments:
        color (int): specific color to use
        count (int): count of metabolites' nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    ranks_summary = prepare_ranks_summary(
        count=count,
        nodes_metabolites=nodes_metabolites
    )
    chart_ranks = plot_three_ranks_parallel_coordinates(
        color=color,
        ranks_summary=ranks_summary,
        fonts=fonts,
        colors=colors
    )
    # Return reference to figure.
    return chart_ranks


def prepare_ranks_summary(
    count=None,
    nodes_metabolites=None
):
    """
    Prepares summary information about ranks of metabolites' nodes.

    arguments:
        count (int): count of metabolites' nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (list<dict<str>): summary of metabolites' nodes' ranks

    """

    # Extract relevant information.
    summary_raw = []
    for node in nodes_metabolites.values():
        identifier = node["identifier"]
        name = node["name"]
        rank_total = node["rank"]
        rank_degree = node["rank_centrality_degree"]
        rank_betweenness = node["rank_centrality_betweenness"]
        record = {
            "identifier": identifier,
            "name": name,
            "rank_total": rank_total,
            "rank_degree": rank_degree,
            "rank_betweenness": rank_betweenness
        }
        summary_raw.append(record)
    # Sort summary by total rank.
    summary_sort = sorted(
        summary_raw,
        key=lambda record: record["rank_total"],
        reverse=False
    )
    # Filter summary by total rank.
    summary_filter = list(filter(
        lambda record: record["rank_total"] < (count + 1),
        summary_sort
    ))
    return summary_filter


def plot_three_ranks_parallel_coordinates(
    color=None,
    ranks_summary=None,
    fonts=None,
    colors=None
):
    """
    Creates a parallel coordinates chart to represent ranks of metabolites'
    nodes.

    arguments:
        color (str): specific color to use
        ranks_summary (list<dict<str>): summary of metabolites' nodes' ranks
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    # Determine minimal and maximal values for all axes.
    minimum, maximum = determine_ranks_summary_extremes(
        ranks_summary=ranks_summary
    )
    # Extract information for each axis.
    categories = ["Degree", "Total", "Betweenness"]
    # Create chart.
    figure = pyplot.figure(
        figsize=(11.811, 15.748),
        tight_layout=True
    )
    axes = pyplot.axes()
    axes.invert_yaxis()
    axes.set_ylabel(
        ylabel="Node Ranks",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["font_one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["font_one"].get_size(),
        labelcolor=colors["black"]
    )
    for category in categories:
        axes.axvline(
            category,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["black"]
        )
    # Create a plot for each value.
    for record in ranks_summary:
        values = [
            record["rank_degree"],
            record["rank_total"],
            record["rank_betweenness"]
        ]
        axes.plot(
            categories,
            values,
            color=colors[color],
            linewidth=2.5
        )
        axes.text(
            categories[1],
            values[1],
            record["name"],
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["font_two"],
            horizontalalignment="center",
            verticalalignment="center"
        )
    # Return reference to figure.
    return figure


def determine_ranks_summary_extremes(ranks_summary=None):
    """
    Determines minimal and maximal values of all ranks in summary.

    arguments:
        ranks_summary (list<dict<str>): summary of metabolites' nodes' ranks

    raises:

    returns:
        (tuple<int, int>): minimal and maximal ranks

    """

    ranks = []
    for record in ranks_summary:
        ranks.append(record["rank_total"])
        ranks.append(record["rank_degree"])
        ranks.append(record["rank_betweenness"])
    minimum = min(ranks)
    maximum = max(ranks)
    return (minimum, maximum)


def plot_names(
    count=None,
    nodes_metabolites=None,
    fonts=None
):
    """
    Creates a chart to represent dominance of metabolites' nodes.

    arguments:
        count (int): count of metabolites' nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        fonts (dict<object>): references to definitions of font properties

    raises:

    returns:
        (object): chart object

    """

    # Determine frequencies of names.
    names_frequencies = {}
    for record in nodes_metabolites.values():
        name = record["name"]
        degree = record["degree"]
        names_frequencies[name] = degree
    # Create word cloud.
    chart_names = wordcloud.WordCloud(
        width=3000,#2000, 4000, 6000,
        height=4000,#1500, 3000, 4500,
        min_font_size=1,
        max_font_size=300,#200, #300, 500
        max_words=count,
        colormap="ocean",#"viridis", "plasma", "ocean", "gist_earth",
        background_color="white",
        prefer_horizontal=0.90,
        relative_scaling=0.75,
        stopwords=[""],
    ).generate_from_frequencies(
        names_frequencies
    )
    # Return reference to figure.
    return chart_names


def plot_volcano(
    records=None,
    threshold_p=None,
    threshold_fold=None,
    label=None,
    fonts=None,
    colors=None
):
    """
    Creates a volcano chart to represent p-values and fold changes.

    arguments:
        records (list<dict): information about values
        threshold_p (float): threshold by p-value
        threshold_fold (float): threshold by fold change
        label (bool): whether to plot labels for points within thresholds
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    # Separate records of values within and without thresholds.
    records_threshold = filter_records_thresholds(
        records=records,
        threshold_p=threshold_p,
        threshold_fold=threshold_fold
    )
    # Create chart.
    figure = pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = pyplot.axes()
    axes.set_xlabel(
        xlabel="log2(fold)",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["font_one"]
    )
    axes.set_ylabel(
        ylabel="-log10(p-value)",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["font_one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["font_one"].get_size(),
        labelcolor=colors["black"]
    )
    # Create lines for thresholds.
    axes.axvline(
        x=math.log(threshold_fold, 2),
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=2.5,
    )
    axes.axvline(
        x=(-1 * math.log(threshold_fold, 2)),
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=2.5,
    )
    axes.axhline(
        y=(-1 * math.log(threshold_p, 10)),
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=2.5,
    )
    # Plot values external to thresholds.
    for record in records_threshold["external"]:
        p_value = record["p_value"]
        p_value_log = (-1 * math.log(p_value, 10))
        fold = record["fold"]
        fold_log = math.log(threshold_fold, 2)
        axes.plot(
            fold_log,
            p_value_log,
            linestyle="",
            marker="o",
            markeredgecolor=colors["blue"],
            markerfacecolor=colors["blue"]
        )
    # Plot values internal to thresholds.
    for record in records_threshold["internal"]:
        name = record["name"]
        p_value = record["p_value"]
        p_value_log = (-1 * math.log(p_value, 10))
        fold = record["fold"]
        fold_log = math.log(threshold_fold, 2)
        axes.plot(
            fold_log,
            p_value_log,
            linestyle="",
            marker="o",
            markeredgecolor=colors["orange"],
            markerfacecolor=colors["orange"]
        )
        axes.text(
            fold_log,
            p_value_log,
            name,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["font_three"],
            horizontalalignment="center",
            verticalalignment="center"
        )
    # Return reference to figure.
    return figure


def filter_records_thresholds(
    records=None,
    threshold_p=None,
    threshold_fold=None
):
    """
    Filters records by whether their values are within or without thresholds.

    arguments:
        records (list<dict): information about values
        threshold_p (float): threshold by p-value
        threshold_fold (float): threshold by fold change

    raises:

    returns:
        (dict<list<dict>>): records within and without thresholds

    """

    # Collect records with values within or without thresholds.
    records_internal = []
    records_external = []
    for record in records:
        p_value = record["p_value"]
        fold = record["fold"]
        match_p = (p_value < threshold_p)
        match_fold = (fold > threshold_fold) or (fold < (1 / threshold_fold))
        if match_p and match_fold:
            records_internal.append(record)
        else:
            records_external.append(record)
    # Compile and return information.
    return {
        "internal": records_internal,
        "external": records_external
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
    path = os.path.join(directory, "plot2")
    utility.confirm_path_directory(path)
    # Degree.
    path_degree_compartments_yes = os.path.join(
        path, "metabolite_degrees_compartments-yes.svg"
    )
    path_degree_compartments_no = os.path.join(
        path, "metabolite_degrees_compartments-no.svg"
    )
    # Rank.
    path_rank_compartments_yes_hubs_yes = os.path.join(
        path, "metabolite_ranks_compartments-yes_hubs-yes.svg"
    )
    path_rank_compartments_yes_hubs_no = os.path.join(
        path, "metabolite_ranks_compartments-yes_hubs-no.svg"
    )
    path_rank_compartments_no_hubs_yes = os.path.join(
        path, "metabolite_ranks_compartments-no_hubs-yes.svg"
    )
    path_rank_compartments_no_hubs_no = os.path.join(
        path, "metabolite_ranks_compartments-no_hubs-no.svg"
    )
    # Name.
    path_name_compartments_yes_hubs_yes = os.path.join(
        path, "metabolite_names_compartments-yes_hubs-yes.png"
    )
    path_name_compartments_yes_hubs_no = os.path.join(
        path, "metabolite_names_compartments-yes_hubs-no.png"
    )
    path_name_compartments_no_hubs_yes = os.path.join(
        path, "metabolite_names_compartments-no_hubs-yes.png"
    )
    path_name_compartments_no_hubs_no = os.path.join(
        path, "metabolite_names_compartments-no_hubs-no.png"
    )
    # Measurement.
    path_measurement_one = os.path.join(path, "measurement_one.svg")
    # Write information to file.
    information["degree_compartments_yes"].savefig(
        path_degree_compartments_yes,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["degree_compartments_no"].savefig(
        path_degree_compartments_no,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["rank_compartments_yes_hubs_yes"].savefig(
        path_rank_compartments_yes_hubs_yes,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["rank_compartments_yes_hubs_no"].savefig(
        path_rank_compartments_yes_hubs_no,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["rank_compartments_no_hubs_yes"].savefig(
        path_rank_compartments_no_hubs_yes,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["rank_compartments_no_hubs_no"].savefig(
        path_rank_compartments_no_hubs_no,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["name_compartments_yes_hubs_yes"].to_file(
        path_name_compartments_yes_hubs_yes
    )
    information["name_compartments_yes_hubs_no"].to_file(
        path_name_compartments_yes_hubs_no
    )
    information["name_compartments_no_hubs_yes"].to_file(
        path_name_compartments_no_hubs_yes
    )
    information["name_compartments_no_hubs_no"].to_file(
        path_name_compartments_no_hubs_no
    )
    information["measurement_one"].savefig(
        path_measurement_one,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
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
    # Define font properties.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()
    # Distributions of nodes' degrees.
    chart_degree_compartments_yes = plot_degrees(
        nodes_hubs_yes=source["nodes_compartments_yes_hubs_yes"],
        nodes_hubs_no=source["nodes_compartments_yes_hubs_no"],
        fonts=fonts,
        colors=colors
    )
    chart_degree_compartments_no = plot_degrees(
        nodes_hubs_yes=source["nodes_compartments_no_hubs_yes"],
        nodes_hubs_no=source["nodes_compartments_no_hubs_no"],
        fonts=fonts,
        colors=colors
    )
    # Ranks of nodes.
    chart_rank_compartments_yes_hubs_yes = plot_ranks(
        color="blue",
        count=10,
        nodes_metabolites=source["nodes_compartments_yes_hubs_yes"],
        fonts=fonts,
        colors=colors
    )
    chart_rank_compartments_yes_hubs_no = plot_ranks(
        color="orange",
        count=10,
        nodes_metabolites=source["nodes_compartments_yes_hubs_no"],
        fonts=fonts,
        colors=colors
    )
    chart_rank_compartments_no_hubs_yes = plot_ranks(
        color="blue",
        count=10,
        nodes_metabolites=source["nodes_compartments_no_hubs_yes"],
        fonts=fonts,
        colors=colors
    )
    chart_rank_compartments_no_hubs_no = plot_ranks(
        color="orange",
        count=10,
        nodes_metabolites=source["nodes_compartments_no_hubs_no"],
        fonts=fonts,
        colors=colors
    )
    # Names of nodes.
    chart_name_compartments_yes_hubs_yes = plot_names(
        count=3000,
        nodes_metabolites=source["nodes_compartments_yes_hubs_yes"],
        fonts=fonts
    )
    chart_name_compartments_yes_hubs_no = plot_names(
        count=3000,
        nodes_metabolites=source["nodes_compartments_yes_hubs_no"],
        fonts=fonts
    )
    chart_name_compartments_no_hubs_yes = plot_names(
        count=3000,
        nodes_metabolites=source["nodes_compartments_no_hubs_yes"],
        fonts=fonts
    )
    chart_name_compartments_no_hubs_no = plot_names(
        count=3000,
        nodes_metabolites=source["nodes_compartments_no_hubs_no"],
        fonts=fonts
    )
    # Volcano plots for metabolomic measurements.
    # source["measurements_one"]
    # include threshold dashed lines for p-value and fold-change...
    # points beyond thresholds should be orange... all others blue...
    chart_measurement_one = plot_volcano(
        records=source["measurements_one"],
        threshold_p=0.01,
        threshold_fold=2.0,
        label=True,
        fonts=fonts,
        colors=colors
    )
    # Compile information.
    information = {
        "degree_compartments_yes": chart_degree_compartments_yes,
        "degree_compartments_no": chart_degree_compartments_no,
        "rank_compartments_yes_hubs_yes": (
            chart_rank_compartments_yes_hubs_yes
        ),
        "rank_compartments_yes_hubs_no": chart_rank_compartments_yes_hubs_no,
        "rank_compartments_no_hubs_yes": chart_rank_compartments_no_hubs_yes,
        "rank_compartments_no_hubs_no": chart_rank_compartments_no_hubs_no,
        "name_compartments_yes_hubs_yes": chart_name_compartments_yes_hubs_yes,
        "name_compartments_yes_hubs_no": chart_name_compartments_yes_hubs_no,
        "name_compartments_no_hubs_yes": chart_name_compartments_no_hubs_yes,
        "name_compartments_no_hubs_no": chart_name_compartments_no_hubs_no,
        "measurement_one": chart_measurement_one,
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
