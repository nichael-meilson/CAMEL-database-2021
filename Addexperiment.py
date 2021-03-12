import urllib

import requests as req
import pandas as pd
import re
import glob
from mutfunc_functionality import extract_files
from mutfunc_functionality import runmutfunc
from mutfunc_functionality import add_column_description
from CellularLocation import *
from Mechismo_functionality import *
import zipfile
path = ""
# Create function to pull all data from excel file, log in and add experiment
# Also includes function to attach mutation data to an experiment

# Only currently works for 1 experiment per entry, can in the future just have it loop through non empty entries
# If there are multiple entries for a field, such as having multiple species in the experiment, these entries need to
# be separated by a comma in the excel template file in order for them to be properly added

# We need to update script so that you can also include a mutation file with it and that file joins the mutation field

def _get_all_mutation_data(mutfile):
    mut_df = pd.read_excel(mutfile, sheet_name='Data', header=4, keep_default_na=False)
    # mut_func_file = runmutfunc(mutfile)
    #TODO remove this later - test purposes only
    mut_func_file = pd.read_csv('test/mutfunc_output_raw.csv')
    # Update our mutation excel file with the Mutfunc results and return it as a new file with appropriately
    # detailed headers
    # updated_mutation_dataframe = extract_files(mut_func_file, mut_df)
    # mechismo_results = run_mechismo(mutfile)
    #TODO remove this later - test purposes only
    mechismo_results = pd.read_csv('test/mechismo_output_raw.csv')
    # Add mechismo results to complete mutation dataframe
    df = pd.merge(mut_func_file, mechismo_results, left_on="Start POS", right_on="Start POS",
                        how='left')
    df = df.fillna('')
    list_of_genes = locations(mutfile)
    # check to see if we can run cell2go with this mutation file, if not we end here
    if list_of_genes == "False":
        return
    # cell2go_results = cello2go(list_of_genes)
    #TODO remove this later - test purposes only
    cell2go_results = pd.read_csv('test/cello_output_raw.csv')
    df = pd.merge(df, cell2go_results, left_on="Start POS", right_on="Start POS", how='left')
    df = df.fillna('')
    # Need to drop duplicates here
    df = df.drop_duplicates()
    return df


def _get_fielddict(val):
    # Create fields dictionary and check for multiple entries by looking for a semicolon
    start = 1
    fielddict = {}
    counter = 1
    while start < len(val)-3:
        if val[start] != '':
            # update for new field of mutation data complete, field ids are tied to excel columns here so need to
            # adjust value to the real id
            # check to see if it is a new field that doesn't go in order, in this case the only one is field 47 which is
            # column 36 in the template and is mutation data complete
            # if start == 36:
            #     fielddict[str(47)] = {}
            #     fielddict[str(47)] = {'new_' + str(counter): val[start]}
            #     start += 1
            #     continue
            fielddict[str(start)] = {}
            # We allow commas in some text fields like comments and major outcomes, even though this doesnt check for
            # actual commas we just let this field go through, might need to update for 35 remarks as well in the future
            if isinstance(val[start], str) and start != 34:
                comma_check = re.split(";", val[start])
                for element in comma_check:
                    fielddict[str(start)]['new_' + str(counter)] = element
                    counter += 1
            else:
                fielddict[str(start)] = {'new_' + str(counter): val[start]}
                counter += 1
        start += 1
    
    return fielddict


def _add_mutations_to_fielddict(fielddict, mut_df):
    # Should start at 36 from fielddict keys
    _start_id = max([int(x) for x in list(fielddict.keys())])
    # Should start at 21 from fielddict values where key = 36
    _counter = int(list(fielddict[str(_start_id)].keys())[0].split('_')[1])

    # Warning: these correspond to field_ids in the SQL db; try not to change them in the DB
    field_mapper = {
        'CHROM_ID': 41,
        'Start POS': 42,
        'End POS': 43,
        'TYPE': 44,
        'REF': 45,
        'ALT': 46,
        'GEN': 47,
        '∆AA': 48,
        'POP': 49,
        'CLON': 50,
        'TIME': 51,
        'FREQ': 52,
        'COM': 53,
        'impact': 54,
        'score': 55,
        'ic': 56,
        'ddg': 57,
        'pdb_id': 58,
        'sequence': 59,
        'accession': 60,
        'modification_type': 61,
        'site_function': 62,
        'function_evidence': 63,
        'predicted_kinase': 64, 
        'probability_loss': 65, 
        'knockout_pvalue': 66, 
        'tf': 67,
        'Category of Mutation': 68,
        'Interactions/Score': 69,
        'Total Interaction Score': 70, 
        'Cello2go probabilities': 71, 
        'Location': 72

    }

    # Parse df to make dicts within dicts within a list
    valid_df = mut_df[field_mapper.keys()]

    mut_list = []
    mut_list.append(fielddict)
    new_counter = _counter+1
    new_start_id = _start_id+1
    for index, values in valid_df.iterrows():
        row_dict = {}
        for i in range(0,len(values)):
            counter_dict = {}
            counter_dict['new_'+str(new_counter)] = values[i]
            row_dict[str(new_start_id)] = counter_dict

            new_counter += 1
            new_start_id += 1

        mut_list.append(row_dict)
        new_counter = _counter+1
        new_start_id = _start_id+1

    # mut_json = str(mut_list)
    # mut_json = "{"+mut_json[1:-1]+"}"

    return mut_list


def _val_serializable(val):
    # Convert data to proper type for updating to database (making everything JSON serializable)
    integer_fields = [3, 16, 17, 26, 27, 28, 29, 32, 33, 38]
    bool_fields = [8, 10, 12, 14, 19, 23, 36]
    double_fields = [30]
    entry = 0
    while entry < len(val):
        if val[entry] != "":
            if entry in integer_fields:
                val[entry] = int(val[entry])
            elif entry in bool_fields:
                val[entry] = bool(val[entry])
            elif entry in double_fields:
                val[entry] = float(val[entry])
        entry += 1
    return val


def get_data_and_add_experiment(file, mutfile =""):

    df = pd.read_excel(file, skiprows=4)
    df = df.fillna("")
    # Take each value that was included as part of the metadata and is not left blank
    val = df.loc[0, :].values.tolist()
    val = _val_serializable(val)
    fielddict = _get_fielddict(val)
    # Get Reference information
    pubmed_id = val[-2]
    if pubmed_id == "":
        pubmed_id = None
    pubmed_url = val[-1]
    # Add a new experiment
    '''
    Adding or updating an element needs a dict that mimics
    the JSON format like a GET request would return

    Field values key/value pairs that do not have a generated ID yet, use a
    random id that is prefixed with 'new_'.
    '''
    # base_url = "https://cameldatabase.com/CAMEL/"
    base_url = "http://localhost:8888/"
    api_url = base_url + "api"
    auth_url = api_url + "/auth/"
    exp_url = api_url + "/experiment"
    attach_url = api_url + '/attachment/'

    # Credentials

    login = 'admin'
    # password = 'oUNr97fbrSr9UVOhj3'
    password = 'password'

    if not password:
        import getpass
        password = getpass.getpass()

    # Get an authentication token
    '''
    All editing operations require a header containing a valid AuthToken.
    A token stays valid for one day.
    '''

    auth_request = req.get(auth_url, auth=(login, password))
    
    token = auth_request.headers['AuthToken']

    # Create the header we are going to send along
    auth_header = {'AuthToken': token}

    new_experiment = {
        'name': df["NAME"][0],  # the only required attribute
        # key value pairs with field id as key
        'fields': fielddict,
        'references': [
            # {
            #     # By default, references need a reference ID and the complete reference data
            #     # with __all reference fields__ (see get results) to do an UPDATE
            #     # This behavior can be changed by the 'action' attribute: 'new' (post, without ref id)
            #     # or 'link' or 'delete' (with existing ref id)
            #     'id': '11954',
            #     'action': 'link'  # link existing paper to this experiment
            # },
            {
                'action': 'new',  # a completely new paper
                'authors': "a list of authors goes here",
                'title': 'this is the title of the new paper',
                'journal': 'Journal Abbr.',
                'year': '2019',
                'pages': '',
                'pubmed_id': pubmed_id,
                'url': pubmed_url
            }
        ]
    }

    # Send the new experiment data
    # It will be added to the database
    # The JSON answer will be the same experiment, but with an assigned ID
    #TODO moved this down
    # answer = req.post(exp_url, headers=auth_header, json=new_experiment).json()
    # exp_id = answer['id']
    # added_exp_url = exp_url + "/" + str(exp_id)

    # we check to see if there is a mutation file attached so we can upload it after the experiment is added, file needs
    # to be .xlsx
    #TODO Fix the 404 later
    # if mutfile != "":
    #     # We upload the file to a temporary location on the server
    #     attachment = {'file': open(mutfile, 'rb')}
    #     resp = req.post(attach_url, files=attachment, headers=auth_header)

    #     # Get the temporary id of the upload
    #     tmp_uuid = resp.json()['uuid']

    #     # Set the attachment field to the tmp id and name the file
    #     # The file will be moved to the correct location
    #     dest_file_name = "Mutation_Data.xlsx"
    #     attach_exp = {
    #         'fields': {
    #             '36': {
    #                 'new_1': {
    #                     'uuid': tmp_uuid,
    #                     'filename': dest_file_name}
    #             }
    #         }
    #     }
    #     resp = req.put(added_exp_url, headers=auth_header, json=attach_exp)
    

    # Last we check to see if a mutation file was attached and then if it qualifies, we have to check the Species
    # to see if it is E. Coli or Yeast(eventually) and then also check to see what the strain is since we only work on
    # the standard most popular strain for both
    # First we have to check the species(starting with just E. Coli)
    if val[1] == "Escherichia coli" and mutfile != "":
        mut_df = _get_all_mutation_data(mutfile)
        mut_df.to_excel("Mutation_results.xlsx", index=False)
        # add_column_description()
        # add this file to the zip file of mutation results
        zip_open = zipfile.ZipFile("Mutation_results_complete.xlsx", 'a')
        zip_open.write("Mutation_results.xlsx")
        zip_open.close()

        # We upload the file to a temporary location on the server
        attachment = {'file': open("Mutation_results_complete.xlsx", 'rb')}
        #TODO Fix gateway issue for attachments
        # resp = req.post(attach_url, files=attachment, headers=auth_header)

        # Add mutation data to experiment post request
        fielddict_mut = _add_mutations_to_fielddict(fielddict, mut_df)

        new_experiment = {
            'name': df["NAME"][0],  # the only required attribute
            # key value pairs with field id as key
            'fields': fielddict_mut,
            'references': [
                # {
                #     # By default, references need a reference ID and the complete reference data
                #     # with __all reference fields__ (see get results) to do an UPDATE
                #     # This behavior can be changed by the 'action' attribute: 'new' (post, without ref id)
                #     # or 'link' or 'delete' (with existing ref id)
                #     'id': '11954',
                #     'action': 'link'  # link existing paper to this experiment
                # },
                {
                    'action': 'new',  # a completely new paper
                    'authors': "a list of authors goes here",
                    'title': 'this is the title of the new paper',
                    'journal': 'Journal Abbr.',
                    'year': '2019',
                    'pages': '',
                    'pubmed_id': pubmed_id,
                    'url': pubmed_url
                }
            ]
        }

        # import pdb; pdb.set_trace()

        answer = req.post(exp_url, headers=auth_header, data=new_experiment).json()
        exp_id = answer['id']
        added_exp_url = exp_url + "/" + str(exp_id)
        # Get the temporary id of the upload
        #TODO change this back
        # tmp_uuid = resp.json()['uuid']
        tmp_uuid = exp_id

        # Set the attachment field to the tmp id and name the file
        # The file will be moved to the correct location
        dest_file_name = "Complete_Mutation_Results.gz"
        attach_exp = {
                'fields': {
                    '44': {
                        'new_1': {
                            'uuid': tmp_uuid,
                            'filename': dest_file_name}
                    }
                }
            }
        resp = req.put(added_exp_url, headers=auth_header, json=attach_exp)

        
        # After we run our script we remove the local version of the files
        # os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results.xlsx")
        # os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results_complete.xlsx")
    else:
        # Field value id's will not be assigned yet, until we request the complete object again
        answer = req.post(exp_url, headers=auth_header, json=new_experiment).json()
        exp_id = answer['id']
        added_exp_url = exp_url + "/" + str(exp_id)
        # answer = req.post(exp_url, headers=auth_header, json=new_experiment).json()
        added_experiment = req.get(added_exp_url).json()


# Function to add mutation data to experiments, needs to be .xlsx

def add_mutation_to_experiment(mutation_file):
    ## URLS
    # base_url = "https://cameldatabase.com/CAMEL/"
    base_url = "http://localhost:8888/"
    api_url = base_url + "api"
    auth_url = api_url + "/auth/"
    exp_url = api_url + "/experiment"
    attach_url = api_url + '/attachment'

    ## Credentials
    login = 'admin'
    # password = 'oUNr97fbrSr9UVOhj3'
    password = 'password'

    if not password:
        import getpass
        password = getpass.getpass()

    ## Get an authentication token
    '''
    All editting operations require a header containing a valid AuthToken.
    A token stays valid for one day.
    '''
    auth_request = req.get(auth_url, auth=(login, password))
    token = auth_request.headers['AuthToken']

    # Create the header we're going to send along
    auth_header = {'AuthToken': token}

    # Get experiment ID from the given file, always needs to be experimentID double _ then file name which should be
    # name of experiment
    file_name = re.split("/", mutation_file)
    file_name = file_name[-1]
    eid = re.split("_", file_name)
    eid = eid[0]

    local_file_name = mutation_file
    exp_id = eid
    added_exp_url = exp_url + '/' + str(exp_id)

    # We upload the file to a temporary location on the server
    attachment = {'file': open(local_file_name, 'rb')}
    resp = req.post(attach_url, files=attachment, headers=auth_header)

    # Get the temporary id of the upload
    tmp_uuid = resp.json()['uuid']

    # Set the attachment field to the tmp id and name the file
    # The file will be moved to the correct location
    dest_file_name = "Mutation_Data.xlsx"
    attach_exp = {
        'fields': {
            '36': {
                'new_1': {
                    'uuid': tmp_uuid,
                    'filename': dest_file_name}
            }
        }
    }
    resp = req.put(added_exp_url, headers=auth_header, json=attach_exp)
    print(resp.json)
    if resp.ok:
        print("Upload successful")
    else:
        print("Upload failed")

    # Now that we attached the mutation data we want to add mutFunc
    # first we have to check the species to see if it is E. coli
    # To do so we use the api and pull the experiment information
    # to check if E. coli or Escherichia coli is a listed species
    experiments = req.get(exp_url + eid).json()
    for key, value in experiments.get("fields").get("1").items():
        if "Escherichia coli" in value:
            # check to see the organism, need to add Yeast or update so its more than just  E. Coli
            mut_df = pd.read_excel(mutation_file, sheet_name='Sheet1', header=4, keep_default_na=False)
            mut_func_file = runmutfunc(mutation_file)
            updated_mutation_dataframe = extract_files(mut_func_file, mut_df)
            mechismo_results = run_mechismo(mutation_file)
            # Add mechismo results to complete mutation dataframe
            df = pd.merge(updated_mutation_dataframe, mechismo_results, left_on="Start POS", right_on="Start POS",
                              how='left')
            df = df.fillna('')
            list_of_genes = locations(mutation_file)
            # check to see if we can run cell2go with this mutation file, if not we end here
            if list_of_genes == "False":
                return
            #TODO comment this out later
            # cell2go_results = run_cello2go(list_of_genes)
            cello2go_results = pd.read_csv('test/cello_output_raw.csv')
            df = pd.merge(df, cell2go_results, left_on="Start POS", right_on="Start POS", how='left')
            df = df.fillna("")
            df.to_excel("Mutation_results.xlsx", index=False)
            add_column_description()
            zip_open = zipfile.ZipFile(mut_func_file, 'a')
            zip_open.write("Mutation_results_complete.xlsx")
            zip_open.close()
            # We upload the file to a temporary location on the server
            attachment = {'file': open(mut_func_file, 'rb')}
            resp = req.post(attach_url, files=attachment, headers=auth_header)

            # Get the temporary id of the upload
            tmp_uuid = resp.json()['uuid']

            # Set the attachment field to the tmp id and name the file
            # The file will be moved to the correct location
            dest_file_name = "Complete_Mutation_Results.gz"
            attach_exp = {
                    'fields': {
                        '44': {
                            'new_1': {
                                'uuid': tmp_uuid,
                                'filename': dest_file_name}
                        }
                    }
                }
            resp = req.put(added_exp_url, headers=auth_header, json=attach_exp)

            # After we run our script we remove the local version of the files
            # os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results.xlsx")
            # os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results_complete.xlsx")


def remove_experiment(eid):

    ## URLS
    # base_url = "https://cameldatabase.com/CAMEL/"
    base_url = "http://localhost:8888/"
    api_url = base_url + "api"
    auth_url = api_url + "/auth/"
    exp_url = api_url + "/experiment"
    attach_url = api_url + '/attachment'

    ## Credentials
    login = 'admin'
    # password = 'oUNr97fbrSr9UVOhj3'
    password = 'password'

    if not password:
        import getpass
        password = getpass.getpass()

    # Get an authentication token
    '''
    All editting operations require a header containing a valid AuthToken.
    A token stays valid for one day.
    '''
    auth_request = req.get(auth_url, auth=(login, password))
    token = auth_request.headers['AuthToken']

    # Create the header we're going to send along
    auth_header = {'AuthToken': token}

    exp_id = eid
    added_exp_url = exp_url + '/' + str(exp_id)
    req.delete(added_exp_url, headers=auth_header)


# remove_experiment(780)
# Have to give file with experiment information and either leave id blank or give a number
# get_data_and_add_experiment('C:/Users/samue/Desktop/Thesis/metadatatemplateUPDATE.xlsx',
#                             "C:/Users/samue/Desktop/Thesis/42C.csv.xlsx")
# get_data_and_add_experiment('C:/Users/samue/Desktop/Thesis/metadatatemplateUPDATE.xlsx')
# Adding experiments from a folder rather than individually
# for fname in glob.glob(path + '\\*'):
#     get_data_and_add_experiment(fname,)
#
# add_mutation_to_experiment('C:/Users/samue/Desktop/Thesis/ALEDB_conversion/MergedExperimentstoUpdate/982_Atsumi S_2010.xlsx')

