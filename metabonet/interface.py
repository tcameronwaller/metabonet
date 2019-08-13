"""
Module to accept parameters from user and execute procedures.

Title:

    metabonet

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
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

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.
import metabonet.metabocurator.clean
import metabonet.metabocurator.collection
import metabonet.metabocurator.conversion
import metabonet.metabocurator.curation
import metabonet.metabocurator.enhancement
import metabonet.metabocurator.extraction
import metabonet.metabocurator.measurement
import metabonet.metabocurator.reconciliation
import metabonet.analysis
import metabonet.candidacy
import metabonet.clean
import metabonet.conversion
import metabonet.measurement
import metabonet.network
import metabonet.plot
import metabonet.utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_interface_parsers():
    """
    Defines and parses arguments from terminal's interface.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define description.
    description = define_general_description()
    # Define epilog.
    epilog = define_general_epilog()
    # Define arguments.
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(title="procedures")
    parser_model = define_model_subparser(subparsers=subparsers)
    parser_network = define_network_subparser(subparsers=subparsers)
    parser_clean = define_clean_subparser(subparsers=subparsers)
    # Parse arguments.
    return parser.parse_args()


def define_general_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Curate model of human metabolism, and define custom networks.

        --------------------------------------------------
    """)
    return description


def define_general_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_model_subparser(subparsers=None):
    """
    Defines subparser for procedures that adapt a model of human metabolism.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_model_description()
    # Define epilog.
    epilog = define_model_epilog()
    # Define parser.
    parser_model = subparsers.add_parser(
        name="model",
        description=description,
        epilog=epilog,
        help="Help for model routine. :)>",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser_model.add_argument(
        "-d", "--directory", dest="directory", type=str, required=True,
        help="Path to root directory for source and product files."
    )
    parser_model.add_argument(
        "-r", "--reconciliation", dest="reconciliation", action="store_true",
        help="Reconcile information from model to MetaNetX."
    )
    parser_model.add_argument(
        "-c", "--collection", dest="collection", action="store_true",
        help="Collect relevant information about metabolic sets and entities."
    )
    parser_model.add_argument(
        "-e", "--extraction", dest="extraction", action="store_true",
        help="Extract information from Human Metabolome Database (HMDB)."
    )
    parser_model.add_argument(
        "-a", "--enhancement", dest="enhancement", action="store_true",
        help="Enhance information about metabolic sets and entities."
    )
    parser_model.add_argument(
        "-u", "--curation", dest="curation", action="store_true",
        help="Curate information about metabolic sets and entities."
    )
    parser_model.add_argument(
        "-v", "--conversion", dest="conversion", action="store_true",
        help="Convert information to formats for export."
    )
    parser_model.add_argument(
        "-m", "--measurement", dest="measurement", action="store_true",
        help="Curate information about measurements of metabolites."
    )
    # Define behavior.
    parser_model.set_defaults(func=evaluate_model_parameters)
    # Return parser.
    return parser_model


def define_model_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        MetaboNet's model procedure

        Curate model of human metabolism.

        --------------------------------------------------
    """)
    return description


def define_model_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        model routine

        --------------------------------------------------
        source information in root directory

        The root directory must contain several sources of information before
        execution of the procedure.

        /root/source/
            recon2m2.xml
            hmdb_metabolites.xml
        /root/source/customization/
            curation_compartments.tsv
            curation_metabolites.tsv
            curation_processes.tsv
            curation_reactions.tsv
            filtration_compartments.tsv
            filtration_processes.tsv
            reconciliation_compartments.tsv
            reconciliation_metabolites.tsv
            simplification_metabolites.tsv
            simplification_reactions.tsv
        /root/source/measurement/metabolomics-workbench_pr000058_st000061/
            analytes.tsv
            measurements.tsv
            samples.tsv
        /root/source/measurement/metabolomics-workbench_pr000305_st000390/
            analytes.tsv
            measurements.tsv
            samples.tsv
        /root/source/measurement/metabolomics-workbench_pr000322_st000412/
            analytes.tsv
            measurements.tsv
            samples.tsv
        /root/source/measurement/metabolomics-workbench_pr000599_st000842/
            analytes.tsv
            measurements.tsv
            samples.tsv

        --------------------------------------------------
        reconciliation

        Reconcile information in human metabolic model to MetaNetX.

        --------------------------------------------------
        collection

        Collect relevant information from metabolic model and from MetaNetX
        about compartments, processes, reactions, and metabolites.

        --------------------------------------------------
        extraction

        Extract relevant information about metabolites from entries in Human
        Metabolome Database (HMDB).

        --------------------------------------------------
        enhancement

        Enhance information about metabolites and reactions.

        --------------------------------------------------
        curation

        Curate information about compartments, processes, reactions, and
        metabolites.

        --------------------------------------------------
        conversion

        Convert information about compartments, processes, reactions, and
        metabolites to versatile formats.

        --------------------------------------------------
        measurement

        Curate information about measurements of metabolites for integration
        with metabolic network.

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_network_subparser(subparsers=None):
    """
    Defines subparser for procedures that define a network.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_network_description()
    # Define epilog.
    epilog = define_network_epilog()
    # Define parser.
    parser_network = subparsers.add_parser(
        name="network",
        description=description,
        epilog=epilog,
        help="Help for network routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser_network.add_argument(
        "-d", "--directory", dest="directory", type=str, required=True,
        help="Path to root directory for source and product files."
    )
    parser_network.add_argument(
        "-y", "--candidacy", dest="candidacy", action="store_true",
        help="Determine candidacy of reactions and metabolites for network."
    )
    parser_network.add_argument(
        "-c", "--compartmentalization", dest="compartmentalization",
        action="store_true", required=False,
        help="Compartmentalize metabolites in candidacy procedure."
    )
    parser_network.add_argument(
        "-s", "--simplification", dest="simplification",
        action="store_true", required=False,
        help=(
            "Simplify specific metabolites and reactions in candidacy " +
            "procedure."
        )
    )
    parser_network.add_argument(
        "-n", "--network", dest="network", action="store_true",
        help=(
            "Define network's nodes and links for candidate reactions and "
            + "metabolites."
        )
    )
    parser_network.add_argument(
        "-p", "--component", dest="component", action="store_true",
        required=False,
        help="Select network's main component in network procedure."
    )
    parser_network.add_argument(
        "-v", "--conversion", dest="conversion", action="store_true",
        help="Convert information to formats for export."
    )
    parser_network.add_argument(
        "-m", "--measurement", dest="measurement", action="store_true",
        help="Integrate information about measurements of metabolites."
    )
    parser_network.add_argument(
        "-a", "--analysis", dest="analysis", action="store_true",
        help="Analyze metabolic network."
    )
    parser_network.add_argument(
        "-t", "--plot", dest="plot", action="store_true",
        help="Generate visual plots from analysis."
    )
    # Define behavior.
    parser_network.set_defaults(func=evaluate_network_parameters)
    # Return parser.
    return parser_network


def define_network_description():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        MetaboNet's network procedure

        Define custom networks from information from human metabolic model.

        --------------------------------------------------
    """)
    return description


def define_network_epilog():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        network routine


        --------------------------------------------------
        candidacy

        Determine candidacy of reactions and metabolites for representation in
        metabolic network.

        Optionally represent compartmental or general metabolites.

        Optionally simplify representations of specific reactions and
        metabolites.

        --------------------------------------------------
        measurement

        Integrate information about measurements of metabolites with candidate
        metabolites for representation in metabolic network.

        --------------------------------------------------
        network

        Define network's nodes and links from reactions and metabolites.

        Optionally select only the network's main component.

        --------------------------------------------------
        conversion

        Convert information about network's elements (nodes and links) to
        formats for transfer to NetworkX and CytoScape.

        --------------------------------------------------
        analysis

        Analyze metabolic network in NetworkX.

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_clean_subparser(subparsers=None):
    """
    Defines subparser for procedures that remove files and directories.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    parser_clean = subparsers.add_parser(
        name="clean",
        description=(
            "Remove all directories and files that these procedures create."
        ),
        help="Help for clean routine."
    )
    parser_clean.add_argument(
        "-d", "--directory", dest="directory", type=str, required=True,
        help="Path to root directory for source and product files."
    )
    # Define behavior.
    parser_clean.set_defaults(func=evaluate_clean_parameters)
    # Return parser.
    return parser_clean


def evaluate_model_parameters(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to model routine ...")
    # Execute procedure.
    if arguments.reconciliation:
        # Report status.
        print("... executing reconciliation procedure ...")
        # Execute procedure.
        metabonet.metabocurator.reconciliation.execute_procedure(
            directory=arguments.directory
        )
    if arguments.collection:
        # Report status.
        print("... executing collection procedure ...")
        # Execute procedure.
        metabonet.metabocurator.collection.execute_procedure(
            directory=arguments.directory
        )
    if arguments.extraction:
        # Report status.
        print("... executing extraction procedure ...")
        # Execute procedure.
        metabonet.metabocurator.extraction.execute_procedure(
            directory=arguments.directory
        )
    if arguments.enhancement:
        # Report status.
        print("... executing enhancement procedure ...")
        # Execute procedure.
        metabonet.metabocurator.enhancement.execute_procedure(
            directory=arguments.directory
        )
    if arguments.curation:
        # Report status.
        print("... executing curation procedure ...")
        # Execute procedure.
        metabonet.metabocurator.curation.execute_procedure(
            directory=arguments.directory
        )
    if arguments.conversion:
        # Report status.
        print("... executing conversion procedure ...")
        # Execute procedure.
        metabonet.metabocurator.conversion.execute_procedure(
            directory=arguments.directory
        )
    if arguments.measurement:
        # Report status.
        print("... executing measurement procedure ...")
        # Execute procedure.
        metabonet.metabocurator.measurement.execute_procedure(
            directory=arguments.directory
        )


def evaluate_network_parameters(arguments):
    """
    Evaluates parameters for network procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to network routine ...")
    # Execute procedure.
    if arguments.candidacy:
        # Report status.
        print("... executing candidacy procedure ...")
        # Execute procedure.
        metabonet.candidacy.execute_procedure(
            compartmentalization=arguments.compartmentalization,
            simplification=arguments.simplification,
            directory=arguments.directory
        )
    if arguments.network:
        # Report status.
        print("... executing network procedure ...")
        # Execute procedure.
        metabonet.network.execute_procedure(
            component=arguments.component,
            directory=arguments.directory
        )
    if arguments.conversion:
        # Report status.
        print("... executing conversion procedure ...")
        # Execute procedure.
        metabonet.conversion.execute_procedure(directory=arguments.directory)
    if arguments.measurement:
        # Report status.
        print("... executing measurement procedure ...")
        # Execute procedure.
        metabonet.measurement.execute_procedure(directory=arguments.directory)
    if arguments.analysis:
        # Report status.
        print("... executing analysis procedure ...")
        # Execute procedure.
        metabonet.analysis.execute_procedure(directory=arguments.directory)
    if arguments.plot:
        # Report status.
        print("... executing plot procedure ...")
        # Execute procedure.
        metabonet.plot.execute_procedure(directory=arguments.directory)


def evaluate_clean_parameters(arguments):
    """
    Evaluates parameters for clean procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    # Report status.
    print("--------------------------------------------------")
    print("... call to clean routine ...")
    # Execute procedure.
    metabonet.metabocurator.clean.execute_procedure(
        directory=arguments.directory
    )
    metabonet.clean.execute_procedure(directory=arguments.directory)



###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # TODO: I want 2 separate procedures: 1. definition, 2. analysis

    # Parse arguments from terminal.
    arguments = define_interface_parsers()
    # Call the appropriate function.
    arguments.func(arguments)


if (__name__ == "__main__"):
    execute_procedure()
