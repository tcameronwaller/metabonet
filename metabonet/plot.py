"""
Author:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 4100, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of MetaboNet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports definition and analysis of custom metabolic networks.
    Copyright (C) 2019 Thomas Cameron Waller

    MetaboNet is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    MetaboNet is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with MetaboNet. If not, see <http://www.gnu.org/licenses/>.
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
import metabonet.utility as utility

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
    nodes_dictionary = {
        "path_true_true": os.path.join(directory, "compartments-true_hubs-true"),
        "path_true_false": os.path.join(directory, "compartments-true_hubs-false"),
        "path_false_true": os.path.join(directory, "compartments-false_hubs-true"),
        "path_false_false": os.path.join(directory, "compartments-false_hubs-false")
        }

    # Nodes
    for key in nodes_dictionary:
        if os.path.exists(nodes_dictionary[key]):

            # Nodes.
            path_nodes_x = os.path.join(
                nodes_dictionary[key], "analysis", "nodes_metabolites.pickle"
            )
            with open(path_nodes_x, "rb") as file_source:
                nodes_x = pickle.load(file_source)
                nodes_dictionary[key] = nodes_x

        else:
            nodes_dictionary[key] = None

    # Measurements.
    measurement_list = []
    path_measurement = os.path.join(directory, "measurement")
    if os.path.exists(path_measurement):

        for file in os.listdir(path_measurement):
            if file.endswith(".pickle"):
                measure_path = os.path.join(path_measurement, file)
                with open(measure_path, "rb") as file_source:
                    measurement_x = pickle.load(file_source)
                    measurement_list.append(measurement_x)
            else:
                pass

    # Check nodes and measurements for non-zero value
    if len(measurement_list) < 1:

        raise Exception("No measurement files found. Exiting...")

    # Compile and return information.
    return {
        "compartmentsTrue_hubsTrue": nodes_dictionary["path_true_true"],
        "compartmentsTrue_hubsFalse": nodes_dictionary["path_true_false"],
        "compartmentsFalse_hubsTrue": nodes_dictionary["path_false_true"],
        "compartmentsFalse_hubsFalse": nodes_dictionary["path_false_false"],
        "measurements": measurement_list
    }


def define_font_properties():
    """
    Defines font properties.

    arguments:

    raises:

    returns:
        (dict<object>): references to definitions of font properties

    """

    # Define font values.
    values_one = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000,
        "weight": 1000,
        "size": 30
    }
    values_two = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 1000,
        "size": 25
    }
    values_three = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 1000,
        "size": 20
    }
    values_four = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 500,
        "size": 15
    }
    # Define font properties.
    properties_one = matplotlib.font_manager.FontProperties(
        family=values_one["family"],
        style=values_one["style"],
        variant=values_one["variant"],
        stretch=values_one["stretch"],
        weight=values_one["weight"],
        size=values_one["size"]
    )
    properties_two = matplotlib.font_manager.FontProperties(
        family=values_two["family"],
        style=values_two["style"],
        variant=values_two["variant"],
        stretch=values_two["stretch"],
        weight=values_two["weight"],
        size=values_two["size"]
    )
    properties_three = matplotlib.font_manager.FontProperties(
        family=values_three["family"],
        style=values_three["style"],
        variant=values_three["variant"],
        stretch=values_three["stretch"],
        weight=values_three["weight"],
        size=values_three["size"]
    )
    properties_four = matplotlib.font_manager.FontProperties(
        family=values_four["family"],
        style=values_four["style"],
        variant=values_four["variant"],
        stretch=values_four["stretch"],
        weight=values_four["weight"],
        size=values_four["size"]
    )
    # Compile and return references.
    return {
        "values": {
            "one": values_one,
            "two": values_two,
            "three": values_three,
            "four": values_four
        },
        "properties": {
            "one": properties_one,
            "two": properties_two,
            "three": properties_three,
            "four": properties_four
        }
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
    nodes_one=None,
    nodes_two=None,
    nodes_three=None,
    nodes_four=None,
    fonts=None,
    colors=None
):
    """
    Creates a chart to represent distributions of degrees of metabolites'
    nodes.

    arguments:
        nodes_one (dict<dict>): information about nodes
        nodes_two (dict<dict>): information about nodes
        nodes_three (dict<dict>): information about nodes
        nodes_four (dict<dict>): information about nodes
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (dict<object>): references to chart objects

    """
    if nodes_one != None and nodes_two != None:
        chart_one_two = plot_degree_distributions(
            nodes_hubs_yes=nodes_one,
            nodes_hubs_no=nodes_two,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_one_two = None

    if nodes_three != None and nodes_four != None:
        chart_three_four = plot_degree_distributions(
            nodes_hubs_yes=nodes_three,
            nodes_hubs_no=nodes_four,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_three_four = None

    # Compile and return information.
    return {
        "one_two": chart_one_two,
        "three_four": chart_three_four
    }


def plot_degree_distributions(
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
    if nodes_hubs_yes == None or nodes_hubs_no == None:
        print("Cannot plot hubs degree comparison.")
        return -1

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
    hist, bin_edges = numpy.histogram(series_one, bins=150)
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
        prop=fonts["properties"]["one"],
        edgecolor=colors["black"]
    )
    axes.set_xlabel(
        xlabel="Node Degree",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel="Count of Nodes",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"]
    )
    return figure


def plot_ranks(
    count=None,
    nodes_one=None,
    nodes_two=None,
    nodes_three=None,
    nodes_four=None,
    fonts=None,
    colors=None
):
    """
    Creates a chart to represent distributions of degrees of metabolites'
    nodes.

    arguments:
        count (int): count of nodes to include in summary
        nodes_one (dict<dict>): information about nodes
        nodes_two (dict<dict>): information about nodes
        nodes_three (dict<dict>): information about nodes
        nodes_four (dict<dict>): information about nodes
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (dict<object>): references to chart objects

    """
    if nodes_one != None:
        chart_one = plot_nodes_ranks(
            color="blue",
            count=count,
            nodes=nodes_one,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_one = None

    if nodes_two != None:
        chart_two = plot_nodes_ranks(
            color="orange",
            count=count,
            nodes=nodes_two,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_two = None

    if nodes_three != None:
        chart_three = plot_nodes_ranks(
            color="blue",
            count=count,
            nodes=nodes_three,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_three = None

    if nodes_four != None:
        chart_four = plot_nodes_ranks(
            color="orange",
            count=count,
            nodes=nodes_four,
            fonts=fonts,
            colors=colors
        )
    else:
        chart_four = None

    # Compile and return information.
    return {
        "one": chart_one,
        "two": chart_two,
        "three": chart_three,
        "four": chart_four
    }


def plot_nodes_ranks(
    color=None,
    count=None,
    nodes=None,
    fonts=None,
    colors=None
):
    """
    Creates a chart to represent ranks of metabolites' nodes.

    arguments:
        color (int): specific color to use
        count (int): count of nodes to include in summary
        nodes (dict<dict>): information about nodes
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): chart object

    """

    ranks_summary = prepare_ranks_summary(
        count=count,
        nodes=nodes
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
    nodes=None
):
    """
    Prepares summary information about ranks of metabolites' nodes.

    arguments:
        count (int): count of metabolites' nodes to include in summary
        nodes (dict<dict>): information about nodes

    raises:

    returns:
        (list<dict<str>): summary of nodes' ranks

    """

    # Extract relevant information.
    summary_raw = []
    for node in nodes.values():
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
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
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
            fontproperties=fonts["properties"]["three"],
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
    nodes_one=None,
    nodes_two=None,
    nodes_three=None,
    nodes_four=None
):
    """
    Creates a chart to represent distributions of degrees of metabolites'
    nodes.

    arguments:
        count (int): count of names to include
        nodes_one (dict<dict>): information about nodes
        nodes_two (dict<dict>): information about nodes
        nodes_three (dict<dict>): information about nodes
        nodes_four (dict<dict>): information about nodes

    raises:

    returns:
        (dict<object>): references to chart objects

    """

    if nodes_one != None:
        print('-> Plotting name cloud for compartments=True and hubs=True...')
        chart_one = plot_names_clouds(
            count=count,
            nodes=nodes_one
        )
    else:
        chart_one = None

    if nodes_two != None:
        print('-> Plotting name cloud for compartments=True and hubs=False...')
        chart_two = plot_names_clouds(
            count=count,
            nodes=nodes_two
        )
    else:
        chart_two = None

    if nodes_three != None:
        print('-> Plotting name cloud for compartments=False and hubs=True...')
        chart_three = plot_names_clouds(
            count=count,
            nodes=nodes_three
        )
    else:
        chart_three = None

    if nodes_four != None:
        print('-> Plotting name cloud for compartments=False and hubs=False...')
        chart_four = plot_names_clouds(
            count=count,
            nodes=nodes_four
        )
    else:
        chart_four = None

    # Compile and return information.
    return {
        "one": chart_one,
        "two": chart_two,
        "three": chart_three,
        "four": chart_four
    }


def plot_names_clouds(
    count=None,
    nodes=None
):
    """
    Creates a chart to represent dominance of metabolites' nodes.

    arguments:
        count (int): count of names to include
        nodes (dict<dict>): information about nodes

    raises:

    returns:
        (object): chart object

    """

    # Determine frequencies of names.
    names_frequencies = {}
    for record in nodes.values():
        name = record["name"]
        # TODO: I can't just use "rank" here because highest rank is 1...
        # TODO: I need a representation of rank that is true to the
        # TODO: distribution of degree and betweenness centrality...
        # TODO: Or I can just use degree for now.
        frequency = record["degree"]
        names_frequencies[name] = frequency
    # Create word cloud.
    chart_names = wordcloud.WordCloud(
        width=3000,#1500, 3000, 4500,
        height=4000,#2000, 4000, 6000,
        min_font_size=1,
        max_font_size=250,#150, 200, 300, 500
        max_words=count,
        colormap="ocean",#"viridis", "plasma", "ocean", "gist_earth",
        background_color="white",
        prefer_horizontal=0.90,
        relative_scaling=1.0,#1.0, 0.75, 0.50
        stopwords=[""],
    ).generate_from_frequencies(
        names_frequencies
    )
    # Return reference to figure.
    return chart_names


def plot_measurements(
    records=None,
    threshold_p=None,
    threshold_fold=None,
    pad=None,
    label=None,
    fonts=None,
    colors=None
):
    """
    Creates charts to represent metabolomic measurements.

    arguments:
        records_one (list<dict>): information about values
        records_two (list<dict>): information about values
        records_three (list<dict>): information about values
        records_four (list<dict>): information about values
        records_five (list<dict>): information about values
        threshold_p (float): threshold by p-value
        threshold_fold (float): threshold by fold change
        label (bool): whether to plot labels for points within thresholds
        pad (float): horizontal pad between point and label
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (dict<object>): references to chart objects

    """
    charts = []

    if len(records) > 0:
        for record_x in records:
            chart_x = plot_volcano(
                records=record_x,
                threshold_p=0.025,
                threshold_fold=1.414,
                label=label,
                pad=pad,
                fonts=fonts,
                colors=colors
            )
            charts.append(chart_x)

        # Compile and return information.
        return charts
    else:
        return None


def plot_volcano(
    records=None,
    threshold_p=None,
    threshold_fold=None,
    label=None,
    pad=None,
    fonts=None,
    colors=None
):
    """
    Creates a volcano chart to represent p-values and fold changes.

    arguments:
        records (list<dict>): information about values
        threshold_p (float): threshold by p-value
        threshold_fold (float): threshold by fold change
        label (bool): whether to plot labels for points within thresholds
        pad (float): horizontal pad between point and label
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
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel="-log10(p-value)",
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
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
        linewidth=2.0,
    )
    axes.axvline(
        x=(-1 * math.log(threshold_fold, 2)),
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=2.0,
    )
    axes.axhline(
        y=(-1 * math.log(threshold_p, 10)),
        xmin=0,
        xmax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=2.0,
    )
    # Plot values external to thresholds.
    for record in records_threshold["external"]:
        p_value = record["p_value"]
        p_value_log = (-1 * math.log(p_value, 10))
        fold = record["fold"]
        fold_log = record["fold_log"]
        axes.plot(
            fold_log,
            p_value_log,
            linestyle="",
            marker="o",
            markersize=15.0,
            markeredgecolor=colors["blue"],
            markerfacecolor=colors["blue"]
        )
    # Plot values internal to thresholds.
    for record in records_threshold["internal"]:
        name = record["name"]
        p_value = record["p_value"]
        p_value_log = (-1 * math.log(p_value, 10))
        fold = record["fold"]
        fold_log = record["fold_log"]
        axes.plot(
            fold_log,
            p_value_log,
            linestyle="",
            marker="o",
            markersize=15.0,
            markeredgecolor=colors["orange"],
            markerfacecolor=colors["orange"]
        )
        if label:
            if fold_log > 0:
                position_horizontal = fold_log - pad
                alignment_horizontal = "right"
            else:
                position_horizontal = fold_log + pad
                alignment_horizontal = "left"
            axes.text(
                position_horizontal,
                (p_value_log),
                name,
                backgroundcolor=colors["white_faint"],
                color=colors["black"],
                fontproperties=fonts["properties"]["four"],
                horizontalalignment=alignment_horizontal,
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
        records (list<dict>): information about values
        threshold_p (float): threshold by p-value
        threshold_fold (float): threshold by fold change

    raises:

    returns:
        (dict<list<dict>>): records within and without thresholds

    """

    # Collect records with values within or without thresholds.
    threshold_fold_log = math.log(threshold_fold, 2)
    records_internal = []
    records_external = []
    for record in records:
        p_value = record["p_value"]
        fold = record["fold"]
        fold_log = record["fold_log"]
        match_p = (p_value < threshold_p)
        match_fold = (
            (fold_log > threshold_fold_log) or
            (fold_log < (-1 * threshold_fold_log))
        )
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
            for all keys except measurements, contains a dictionary of several plots or None if no plot available
            for measurements, list of plots or None

    raises:

    returns:

    """

    # Specify directories and files.
    path = os.path.join(directory, "plot")
    utility.confirm_path_directory(path)

    # Degree.
    path_degree_one_two = os.path.join(path, "metabolite_degrees_compartments-true.svg")
    path_degree_three_four = os.path.join(path, "metabolite_degrees_compartments-false.svg")

    # Rank.
    path_rank_one = os.path.join(path, "metabolite_ranks_compartments-true_hubs-true.svg")
    path_rank_two = os.path.join(path, "metabolite_ranks_compartments-true_hubs-false.svg")
    path_rank_three = os.path.join(path, "metabolite_ranks_compartments-false_hubs-true.svg")
    path_rank_four = os.path.join(path, "metabolite_ranks_compartments-false_hubs-false.svg")

    # Name.
    path_name_one = os.path.join(path, "metabolite_names_compartments-true_hubs-true.png")
    path_name_two = os.path.join(path, "metabolite_names_compartments-true_hubs-false.png")
    path_name_three = os.path.join(path, "metabolite_names_compartments-false_hubs-true.png")
    path_name_four = os.path.join(path, "metabolite_names_compartments-false_hubs-false.png")

    # Write information to file.
    if information["charts_degrees"]["one_two"] != None:
        information["charts_degrees"]["one_two"].savefig(
            path_degree_one_two,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )
    if information["charts_degrees"]["three_four"] != None:
        information["charts_degrees"]["three_four"].savefig(
            path_degree_three_four,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )

    if information["charts_ranks"]["one"] != None:
        information["charts_ranks"]["one"].savefig(
            path_rank_one,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )
    if information["charts_ranks"]["two"] != None:
        information["charts_ranks"]["two"].savefig(
            path_rank_two,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )
    if information["charts_ranks"]["three"] != None:
        information["charts_ranks"]["three"].savefig(
            path_rank_three,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )
    if information["charts_ranks"]["four"] != None:
        information["charts_ranks"]["four"].savefig(
            path_rank_four,
            format="svg",
            #dpi=600,
            facecolor="w",
            edgecolor="w",
            transparent=False
        )

    if information["charts_names"]["one"] != None:
        information["charts_names"]["one"].to_file(path_name_one)

    if information["charts_names"]["two"] != None:
        information["charts_names"]["two"].to_file(path_name_two)

    if information["charts_names"]["three"] != None:
        information["charts_names"]["three"].to_file(path_name_three)

    if information["charts_names"]["four"] != None:
        information["charts_names"]["four"].to_file(path_name_four)

    if len(information["charts_measurements"]) > 0:

        figure_number = 0
        figure_total = len(information["charts_measurements"])
        counter = 1

        for x in range(figure_total):
            if information["charts_measurements"][x] != None:
                path_measurement_x = os.path.join(path, "measurements_" + str(counter) + ".svg")
                information["charts_measurements"][x].savefig(
                    path_measurement_x,
                    format="svg",
                    #dpi=600,
                    facecolor="w",
                    edgecolor="w",
                    transparent=False
                )

            figure_number += 1
            counter += 1

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
    print('Plotting node degree plot(s)...')
    charts_degrees = plot_degrees(
        nodes_one=source["compartmentsTrue_hubsTrue"],
        nodes_two=source["compartmentsTrue_hubsFalse"],
        nodes_three=source["compartmentsFalse_hubsTrue"],
        nodes_four=source["compartmentsFalse_hubsFalse"],
        fonts=fonts,
        colors=colors
    )
    # Ranks of nodes.
    print('Plotting rank degree plot(s)...')
    charts_ranks = plot_ranks(
        count=10,
        nodes_one=source["compartmentsTrue_hubsTrue"],
        nodes_two=source["compartmentsTrue_hubsFalse"],
        nodes_three=source["compartmentsFalse_hubsTrue"],
        nodes_four=source["compartmentsFalse_hubsFalse"],
        fonts=fonts,
        colors=colors
    )
    # Names of nodes.
    print('Plotting node names plot(s)...')
    charts_names = plot_names(
        count=5000,
        nodes_one=source["compartmentsTrue_hubsTrue"],
        nodes_two=source["compartmentsTrue_hubsFalse"],
        nodes_three=source["compartmentsFalse_hubsTrue"],
        nodes_four=source["compartmentsFalse_hubsFalse"],
    )
    # Metabolomic measurements.
    print('Plotting metabolic measurement plot(s)...')
    charts_measurements = plot_measurements(
        records=source["measurements"],
        threshold_p=0.01, # 0.05, 0.01
        threshold_fold=1.414, # 1.2, 2.0
        label=True,
        pad=0.1,
        fonts=fonts,
        colors=colors
    )
    # Compile information.
    information = {
        "charts_degrees": charts_degrees, #Dict of network reps with plot or None
        "charts_ranks": charts_ranks, #Dict of network reps with plot or None
        "charts_names": charts_names, #Dict of network reps with plot or None
        "charts_measurements": charts_measurements #List or None
    }
    #Write product information to file.
    print('Writing plots...')
    write_product(directory=directory, information=information)
