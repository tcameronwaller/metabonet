"""
Module to accept parameters from user and execute procedures.

Title:

    metabocurator

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

import reconciliation
import extraction
import enhancement

#dir()
#importlib.reload()

###############################################################################
# Functionality

def parse_arguments():
    """
    Defines and interprets arguments from terminal.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define arguments.
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            --------------------------------------------------

            Curate model of metabolism for definition of networks.

            --------------------------------------------------
        """),
        epilog=textwrap.dedent("""\

            --------------------------------------------------

            reconciliation procedure

            source
            1. model of human metabolism in
               Systems Biology Markup Language (SBML) format
            2. curation information about metabolites in
               text table with tab delimiters

            product
            ...

            --------------------------------------------------

            adaptation procedure

            source
            ...

            product
            ...

            --------------------------------------------------
        """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o", "--origin", dest="origin", type=str, required=True,
        help="Directory of source files."
    )
    parser.add_argument(
        "-d", "--destination", dest="destination", type=str, required=True,
        help="Directory for product files."
    )
    procedure = parser.add_mutually_exclusive_group(required=True)
    procedure.add_argument(
        "-r", "--reconciliation", dest="reconciliation", action="store_true",
        help="Reconcile model to MetaNetX."
    )
    procedure.add_argument(
        "-a", "--adaptation", dest="adaptation", action="store_true",
        help="Adapt model from MetaNetX."
    )
    parser.add_argument(
        "-x", "--clean", dest="clean", action="store_true",
        help="Clean intermediate files."
    )
    # Parse arguments.
    return parser.parse_args()


def evaluate_source(directory=None):
    """
    Evaluates source files.

    arguments:

        directory (object): arguments from terminal

    raises:

    returns:
        (bool): whether arguments are adequate

    """

    # TODO: Make sure necessary input files are available, etc...
    return True

###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    raises:

    returns:

    """

    # Parse arguments from terminal.
    arguments = parse_arguments()
    # Execute procedure.
    if arguments.reconciliation:
        # Report status.
        print("... executing reconciliation procedure ...")
        # Execute reconciliation procedure.
        reconciliation.execute_procedure(
            origin=arguments.origin,
            destination=arguments.destination,
            clean=arguments.clean
        )
    elif arguments.adaptation:
        # Report status.
        print("... executing adaptation procedure ...")
        # Execute adaptation procedure.
        extraction.execute_procedure(
            origin=arguments.origin,
            destination=arguments.destination,
            clean=arguments.clean
        )
        if False:
            enhancement.execute_procedure(
                origin=arguments.origin,
                destination=arguments.destination,
                clean=arguments.clean
            )
            # TODO: spit out a report after enhancement for sake of curation...


            curation.execute_procedure(
                origin=arguments.origin,
                destination=arguments.destination,
                clean=arguments.clean
            )
            conversion.execute_procedure(
                origin=arguments.origin,
                destination=arguments.destination,
                clean=arguments.clean
            )


if (__name__ == "__main__"):
    execute_procedure()
