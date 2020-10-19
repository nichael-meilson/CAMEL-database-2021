import urllib

import requests as req
import pandas as pd
import re
import glob
from MutFunc_functionality import extract_files
from MutFunc_functionality import runmutfunc
from MutFunc_functionality import add_column_description
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


def get_data_and_add_experiment(file, mutfile =""):

    df = pd.read_excel(file, skiprows=4)
    df = df.fillna("")
    # Take each value that was included as part of the metadata and is not left blank
    val = df.loc[0, :].values.tolist()
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
            if start == 36:
                fielddict[str(47)] = {}
                fielddict[str(47)] = {'new_' + str(counter): val[start]}
                start += 1
                continue
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
    print(fielddict)
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
    base_url = "https://cameldatabase.com/CAMEL/"
    api_url = base_url + "api"
    auth_url = base_url + "auth"
    exp_url = api_url + "/experiment"
    field_url = api_url + "/field"
    ref_url = api_url + "/reference"
    attach_url = api_url + '/attachment'
    # Credentials

    login = ''
    password = ''

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
    answer = req.post(exp_url, headers=auth_header, json=new_experiment).json()
    exp_id = answer['id']
    added_exp_url = exp_url + "/" + str(exp_id)

    # we check to see if there is a mutation file attached so we can upload it after the experiment is added, file needs
    # to be .xlsx
    if mutfile != "":
        # We upload the file to a temporary location on the server
        attachment = {'file': open(mutfile, 'rb')}
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
    # Field value id's will not be assigned yet, until we request the complete object again
    added_experiment = req.get(added_exp_url).json()

    # Last we check to see if a mutation file was attached and then if it qualifies, we have to check the Species
    # to see if it is E. Coli or Yeast(eventually) and then also check to see what the strain is since we only work on
    # the standard most popular strain for both
    # First we have to check the species(starting with just E. Coli)
    if val[1] == "Escherichia coli" and mutfile != "":
        mut_df = pd.read_excel(mutfile, sheet_name='Sheet1', header=4, keep_default_na=False)
        mut_func_file = runmutfunc(mutfile)
        # Update our mutation excel file with the Mutfunc results and return it as a new file with appropriately
        # detailed headers
        updated_mutation_dataframe = extract_files(mut_func_file, mut_df)
        mechismo_results = run_mechismo(mutfile)
        # Add mechismo results to complete mutation dataframe
        df = pd.merge(updated_mutation_dataframe, mechismo_results, left_on="Start POS", right_on="Start POS",
                          how='left')
        df = df.fillna('')
        list_of_genes = locations(mutfile)
        # check to see if we can run cell2go with this mutation file, if not we end here
        if list_of_genes == "False":
            return
        cell2go_results = cello2go(list_of_genes)
        df = pd.merge(df, cell2go_results, left_on="Start POS", right_on="Start POS", how='left')
        df = df.fillna('')
        # Need to drop duplicates here
        df = df.drop_duplicates()
        df.to_excel("Mutation_results.xlsx", index=False)
        add_column_description()
        # add this file to the zip file of mutation results
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
        os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results.xlsx")
        os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results_complete.xlsx")


# Function to add mutation data to experiments, needs to be .xlsx

def add_mutation_to_experiment(mutation_file):
    ## URLS
    base_url = "https://cameldatabase.com/CAMEL/"
    api_url = base_url + "api"
    auth_url = base_url + "auth"

    exp_url = api_url + "/experiment"
    attach_url = api_url + '/attachment'

    ## Credentials
    login = ''
    password = ''

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
    if resp.ok:
        print("Upload successful")
    else:
        print("Upload failed")

    # Now that we attached the mutation data we want to add mutFunc
    # first we have to check the species to see if it is E. coli
    # To do so we use the api and pull the experiment information
    # to check if E. coli or Escherichia coli is a listed species
    experiments = req.get("https://cameldatabase.com/CAMEL/api/experiment/" + eid).json()
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
            cell2go_results = cello2go(list_of_genes)
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
            os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results.xlsx")
            os.remove("C:\\Users\\samue\\PycharmProjects\\Thesis\\Mutation_results_complete.xlsx")


def remove_experiment(eid):

    # URLS
    base_url = "https://cameldatabase.com/CAMEL/"
    api_url = base_url + "api"
    auth_url = base_url + "auth"

    exp_url = api_url + "/experiment"
    field_url = api_url + "/field"
    ref_url = api_url + "/reference"

    # Credentials
    login = ''
    password = ''

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
add_mutation_to_experiment('C:/Users/samue/Desktop/Thesis/ALEDB_conversion/MergedExperimentstoUpdate/982_Atsumi S_2010.xlsx')

