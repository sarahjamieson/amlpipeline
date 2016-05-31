import pandas as pd


class ParseSampleSheet(object):
    """Parses a sample sheet (CSV format) into two Python dictionaries, one for header details and one for sample details.

        :param csv_file: sample sheet

       Notes:
           Functions in this class are based on a standard sample sheet layout with the following header attributes:
            - Header
            - Manifests
            - Reads
            - Settings
            - Data
    """
    def __init__(self, csv_file):
        self.csv = csv_file

    def parse_sample_sheet(self):
        # -----------------------------------------------------------------------------------------------------------
        # 1) Set some variables so we can use outside loops
        # -----------------------------------------------------------------------------------------------------------
        header_index = 0
        data_index = 0
        manifest_index = 0
        read_index = 0
        settings_index = 0
        df_run_data_temp = pd.DataFrame([])
        df_run_data_final = pd.DataFrame(columns=['Property', 'Value'])  # this allows for easy appending later
        # -----------------------------------------------------------------------------------------------------------
        # 2) Parse sample sheet into pandas dataframe
        # -----------------------------------------------------------------------------------------------------------
        df_sample_sheet = pd.read_csv(self.csv, header=None)
        # -----------------------------------------------------------------------------------------------------------
        # 3) Get indexes where these details are
        for column in df_sample_sheet:
            for row_index, row in df_sample_sheet.iterrows():
                if row[column] == '[Data]':
                    data_index = row_index
                    df_headings = df_sample_sheet.ix[:data_index - 2, 0:1]  # Put all header info into a separate dataframe
                    df_headings.columns = ['Property', 'Value']
                elif row[column] == '[Header]':
                    header_index = row_index
                elif row[column] == '[Manifests]':
                    manifest_index = row_index
                elif row[column] == '[Reads]':
                    read_index = row_index
                elif row[column] == '[Settings]':
                    settings_index = row_index
                else:
                    pass
        # ----------------------------------------------------------------------------------------------------------
        # 4) Look at header info first: separate the header types and modify to correctly re-merge later.
        # ----------------------------------------------------------------------------------------------------------
        # [Header]
        df_headers = df_run_data_temp.ix[header_index + 1:manifest_index - 1]
        # [Manifests]
        df_manifests = df_run_data_temp.ix[manifest_index + 1:read_index - 2]
        for row_index, row in df_manifests.iterrows():
            row['Property'] = 'Manifest ' + row['Property']
        # [Reads]
        df_reads = df_run_data_temp.ix[read_index + 1:settings_index - 2]
        read_list = []
        for row_index, row in df_reads.iterrows():
            read_list.append(row['Property'])
        # [Settings]
        df_settings = df_run_data_temp.ix[settings_index + 1:]
        # Combine all
        df_run_data_final = df_run_data_final.append(df_headers)
        df_run_data_final = df_run_data_final.append(df_manifests)
        df_run_data_final = df_run_data_final.append({'Property': 'Reads', 'Value': read_list}, ignore_index=True)
        df_run_data_final = df_run_data_final.append(df_settings)
        df_run_data_final = df_run_data_final.reset_index(drop=True)
        # Convert to dictionary, set_index avoids the index being used a key
        run_dict = df_run_data_final.set_index('Property')['Value'].to_dict()
        # ----------------------------------------------------------------------------------------------------------
        # 5) Now look at sample data: extract lab numbers and transpose dataframe to make dictionary work per patient.
        # ----------------------------------------------------------------------------------------------------------
        df_data = df_sample_sheet.ix[data_index + 1:]
        df_data = df_data.reset_index(drop=True)
        # Change column names
        df_data.columns = df_data.iloc[0]
        df_data = df_data.reindex(df_data.index.drop(0))
        # Drop any columns with "NaN" all the way through
        df_data = df_data.dropna(axis=1, how='all')
        # Use lab numbers as column headings and initial key in dictionary
        sample_id_list = []
        for row_index, row in df_data.iterrows():
            sample_id_list.append(row['Sample_Name'][3:12])
        df_data_trans = df_data.transpose()
        df_data_trans.columns = sample_id_list
        # Convert to dictionary
        sample_dict = df_data_trans.to_dict()
        # ----------------------------------------------------------------------------------------------------------
        return run_dict, sample_dict
