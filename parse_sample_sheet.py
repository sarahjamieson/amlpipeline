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
                    df_run_data_temp = df_sample_sheet.ix[:data_index - 2, 0:1]  # Put all header info into a separate df
                    df_run_data_temp.columns = ['Property', 'Value']
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

    '''
    Method 2:

    def get_run_info(csv_file):
    iem = ''
    investigator = ''
    experiment = ''
    run_date = ''
    workflow = ''
    app = ''
    assay = ''
    description = ''
    chemistry = ''
    worksheet = ''
    manifest = ''
    reads = ''
    data_index = 0
    read_index = 0
    manifest_index = 0
    read1 = 0
    read2 = 0
    sample_dict = {}

    with open(csv_file, 'r') as c:
        reader = csv.reader(c, delimiter=',')
        for i, row in enumerate(reader):
            if row[0] == 'IEMFileVersion':
                iem = row[1]
            elif row[0] == "Investigator Name":
                investigator = row[1]
            elif row[0] == 'Experiment Name':
                experiment = row[1]
            elif row[0] == 'Date':
                run_date = row[1]
            elif row[0] == 'Workflow':
                workflow = row[1]
            elif row[0] == 'Application':
                app = row[1]
            elif row[0] == 'Assay':
                assay = row[1]
            elif row[0] == 'Description':
                description = row[1]
            elif row[0] == 'Chemistry':
                chemistry = row[1]
            elif row[0] == 'worksheet':
                worksheet = row[1]
            elif row[0] == '[Manifests]':
                manifest_index = i
            elif row[0] == '[Reads]':
                read_index = i
            elif row[0] == '[Data]':
                data_index = i
            else:
                pass
            if i == (read_index + 1):
                read1 = row[0]
            if i == (read_index + 2):
                read2 = row[0]
            if i == (manifest_index + 1):
                manifest = row[1]
            reads = "(%s,%s)" % (read1, read2)

    run_dict = {
        "IEM": iem,
        "Investigator": investigator,
        "Experiment": experiment,
        "Date": run_date,
        "Workflow": workflow,
        "Application": app,
        "Assay": assay,
        "Description": description,
        "Chemistry": chemistry,
        "worksheet": worksheet,
        "Manifest": manifest,
        "Reads": reads
    }

    df_sample_sheet = pd.read_csv(csv_file, header=None)
    df_data = df_sample_sheet.ix[data_index + 1:]
    for row_index, row in df_data.iterrows():
        lab_id = str(row[1])[3:12]
        sample_id = row[0]
        name = row[1]
        plate = row[2]
        well = row[3]
        index1 = row[5]
        index2 = row[7]
        sample_manifest = row[8]
        project = row[10]

        sample_dict[lab_id] = {
            "Sample_id": sample_id,
            "Name": name,
            "Plate": plate,
            "Well": well,
            "Index1": index1,
            "Index2": index2,
            "Manifest": sample_manifest,
            "Project": project
        }

    print run_dict, sample_dict
    '''
