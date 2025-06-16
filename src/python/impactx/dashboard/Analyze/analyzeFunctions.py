"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import pandas as pd

from .. import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class AnalyzeFunctions:
    """
    Helper functions for
    preparing the contents for the 'Analyze' page
    """

    # -----------------------------------------------------------------------------
    # Functions for Beam Characteristic and Ref particle data table
    # -----------------------------------------------------------------------------

    @staticmethod
    def load_data(file_path):
        """
        Reads data from the provided file path.
        :param file_path: The path to the file to be read.
        :return: A DataFrame containing the data from the file.
        """

        df = pd.read_csv(file_path, sep=" ")
        return df

    @staticmethod
    def convert_to_dict(combined_data):
        """
        Converts data to dictionary format
        for Vuetify data table.
        :param combined_data: The DataFrame to be converted.
        :return: A tuple containing the data as a list of dictionaries and the headers.
        """

        dictionary = combined_data.to_dict(orient="records")
        columns = combined_data.columns
        headers = [
            {"title": column.strip(), "key": column.strip()} for column in columns
        ]
        return dictionary, headers

    @staticmethod
    def combine_files(file1_name, file2_name):
        """
        Merges two files together.
        :param file1_name: The name of the first file.
        :param file2_name: The name of the second file.
        :return: A DataFrame containing the merged data from the two files.
        """

        file1 = AnalyzeFunctions.load_data(file1_name)
        file2 = AnalyzeFunctions.load_data(file2_name)
        return pd.merge(file1, file2, how="outer")

    @staticmethod
    def filter_headers() -> list:
        """
        Retrieves only user-selected headers.
        
        :return: A list of filtered headers based on user selection.
        """

        all_headers = state.all_headers
        selected_headers = state.selected_headers

        filtered_headers = []
        for selectedHeader in all_headers:
            if selectedHeader["key"] in selected_headers:
                filtered_headers.append(selectedHeader)
        return filtered_headers

    @staticmethod
    def filter_data() -> list:
        """
        Retrieves data for user-selected headers.

        :return: A list of filtered data based on user selection.
        """

        all_data = state.all_data
        selected_headers = state.selected_headers
 
        filtered_data = []
        for row in all_data:
            filtered_row = {}
            for key, value in row.items():
                if key in selected_headers:
                    filtered_row[key] = value
            filtered_data.append(filtered_row)
        return filtered_data
