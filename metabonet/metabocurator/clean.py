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

# Relevant

# Custom
import metabonet.utility as utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to convert information about metabolic
    entities and sets to versatile formats.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Specify directories and files.
    path_reconciliation = os.path.join(directory, "reconciliation")
    path_collection = os.path.join(directory, "collection")
    path_extraction = os.path.join(directory, "extraction")
    path_enhancement = os.path.join(directory, "enhancement")
    path_curation = os.path.join(directory, "curation")
    path_conversion = os.path.join(directory, "conversion")
    path_measurement = os.path.join(directory, "measurement")
    # Remove directories and files.
    # Reconciliation.
    if False:
        utility.remove_file(path_reconciliation, "recon2m2_reconciliation.xml")
    # Collection.
    utility.remove_file(os.path.join(path_collection, "compartments.pickle"))
    utility.remove_file(os.path.join(path_collection, "processes.pickle"))
    utility.remove_file(os.path.join(path_collection, "reactions.pickle"))
    utility.remove_file(os.path.join(path_collection, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_collection, "reactions.tsv"))
    utility.remove_file(os.path.join(path_collection, "metabolites.tsv"))
    utility.remove_empty_directory(path_collection)
    # Extraction.
    utility.remove_file(os.path.join(path_extraction, "hmdb_summary.pickle"))
    utility.remove_file(os.path.join(path_extraction, "hmdb_summary.tsv"))
    utility.remove_empty_directory(path_extraction)
    # Enhancement.
    utility.remove_file(os.path.join(path_enhancement, "compartments.pickle"))
    utility.remove_file(os.path.join(path_enhancement, "processes.pickle"))
    utility.remove_file(os.path.join(path_enhancement, "reactions.pickle"))
    utility.remove_file(os.path.join(path_enhancement, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_enhancement, "reactions.tsv"))
    utility.remove_file(os.path.join(path_enhancement, "metabolites.tsv"))
    utility.remove_file(os.path.join(path_enhancement, "reactions_filter.tsv"))
    utility.remove_empty_directory(path_enhancement)
    # Curation.
    utility.remove_file(os.path.join(path_curation, "compartments.pickle"))
    utility.remove_file(os.path.join(path_curation, "processes.pickle"))
    utility.remove_file(os.path.join(path_curation, "reactions.pickle"))
    utility.remove_file(os.path.join(path_curation, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_curation, "reactions.tsv"))
    utility.remove_file(os.path.join(path_curation, "metabolites.tsv"))
    utility.remove_empty_directory(path_curation)
    # Conversion.
    utility.remove_file(os.path.join(path_conversion, "compartments.pickle"))
    utility.remove_file(os.path.join(path_conversion, "processes.pickle"))
    utility.remove_file(os.path.join(path_conversion, "reactions.pickle"))
    utility.remove_file(os.path.join(path_conversion, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_conversion, "compartments.tsv"))
    utility.remove_file(os.path.join(path_conversion, "processes.tsv"))
    utility.remove_file(os.path.join(path_conversion, "reactions.tsv"))
    utility.remove_file(os.path.join(path_conversion, "metabolites.tsv"))
    utility.remove_file(os.path.join(path_conversion, "dymetabonet.json"))
    utility.remove_empty_directory(path_conversion)
    # Measurement.
    utility.remove_file(os.path.join(path_measurement, "study_one.pickle"))
    utility.remove_file(os.path.join(path_measurement, "study_two.pickle"))
    utility.remove_file(os.path.join(path_measurement, "study_three.pickle"))
    utility.remove_file(os.path.join(path_measurement, "study_four.pickle"))
    utility.remove_file(os.path.join(path_measurement, "study_five.pickle"))
    utility.remove_file(os.path.join(path_measurement, "study_one.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_two.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_three.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_four.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_five.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_one_report.tsv"))
    utility.remove_file(os.path.join(path_measurement, "study_two_report.tsv"))
    utility.remove_file(os.path.join(
        path_measurement, "study_three_report.tsv"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_four_report.tsv"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_five_report.tsv"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_one_metaboanalyst.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_two_metaboanalyst.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_three_metaboanalyst.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_four_metaboanalyst.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_five_metaboanalyst.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_one_metaboanalyst_pair.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_two_metaboanalyst_pair.txt"
    ))
    utility.remove_file(os.path.join(
        path_measurement, "study_five_metaboanalyst_pair.txt"
    ))
    utility.remove_empty_directory(path_measurement)
