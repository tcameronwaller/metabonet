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

# Standard
import argparse
import textwrap

# Relevant

# Custom
import candidacy
import network

#dir()
#importlib.reload()

###############################################################################
# Functionality

def define_parse_arguments():
    """
    Defines and parses arguments from terminal.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define arguments.
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            --------------------------------------------------

            Define custom networks from metabolic models.

            --------------------------------------------------
        """),
        epilog=textwrap.dedent("""\

            --------------------------------------------------

            --------------------------------------------------
        """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-d", "--directory", dest="directory", type=str, required=True,
        help="Path to directory for source and product files."
    )
    parser.add_argument(
        "-c", "--compartmentalize", dest="compartmentalization",
        action="store_true", required=False,
        help="Compartmentalize metabolites."
    )
    parser.add_argument(
        "-x", "--clean", dest="clean", action="store_true", required=False,
        help="Clean intermediate files."
    )
    # Parse arguments.
    return parser.parse_args()


def evaluate_source(directory=None):
    """
    Evaluates source files.

    arguments:

        directory (object): arguments from terminal

    returns:
        (bool): whether arguments are adequate

    raises:

    """


    # TODO: Make sure necessary input files are available, etc...
    return True

###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # Parse arguments from terminal.
    arguments = define_parse_arguments()
    # Execute procedure.
    candidacy.execute_procedure(
        compartmentalization=arguments.compartmentalization,
        directory=arguments.directory
    )
    network.execute_procedure(directory=arguments.directory)
    # TODO: maybe include a separate module for measurements...
    if arguments.clean:
        # Report status.
        print("... executing clean procedure ...")
        # Execute clean procedure.
        clean.execute_procedure(
            directory=arguments.directory
        )


if (__name__ == "__main__"):
    execute_procedure()
