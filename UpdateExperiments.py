import requests as req
import pandas as pd
# Thomas' script as a guide to refer to how the api works
# This script will update experiments by saving the experiment number it originated as, deleting the experiment
# and re-adding the experiment with the old experiment number
# Not the current goal to use this script and to instead use addexperiment for updating and adding
'''
This is a simple script to demonstrate the editing API 
of CAMEL. If you run it one line at a time, you can see
how different entities are added, updated and removed.
At the end of the script, the db should be back in its original
state.
Request types:
'GET': retrieve one or more entities as described in the README
'POST': inserts a new experiment or field
'PUT': update an existing experiment or field given the id
'DELETE': remove an existing experiment, field or reference given the id
References can be added and updated as a part of an experiment.
Removing a field will automatically remove all data of that type, just like
removing an experiment will also delete its field data. Orphan references
(references that are not linked to any experiment anymore) are
automatically removed as well.
'''
# first we pull data from our excel file
df = pd.read_excel('C:/Users/samue/Desktop/Thesis/metadatatemplate.xlsx', skiprows=4)
df = df.fillna("")
# Take each value that was included as part of the metadata and is not left blank
val = df.loc[0, :].values.tolist()

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
All editing operations require a header containing a valid AuthToken.
A token stays valid for one day.
'''
auth_request = req.get(auth_url, auth=(login, password))
token = auth_request.headers['AuthToken']

# Create the header we are going to send along
auth_header = {'AuthToken': token}

#WE need to know the id of the experiment that we want to update
# This should be given by the people who provide the information that needs to be updated or else we need to search for it in the
# database
exp_id = ""

added_exp__url = exp_url + "/" + str(exp_id)

# Field value id's will not be assigned yet, until we request the complete object again
added_experiment = req.get(added_exp__url).json()

# Lets update the Species field. We only need the changing fields and the ids of their values.
species_value_ids = list(added_experiment['fields']['1'].keys())

# Create fields dictionary
start = 1
update_dict = {}
counter = 1
while start < len(val)-3:
    if val[start] != '':
        update_dict[str(start)] = {added_experiment['fields'][start].keys(): val[start]}
        counter += 1
    start += 1

# Create reference dictionary depending on if/what reference parts need to be updated.
# Can make a loop if fields get larger
refdict = {}
if val[-2] != "":
    refdict['pubmed_id'] = val[-2]
if val[-1] != "":
    refdict['url'] = val[-1]

# Check to see if reference section needs to be updated as well,
# needs to change if we add more reference columns to metadata
if val[-1] != "" or val[-2]:
    update_experiment = {
        'fields': update_dict,
        'references': [refdict]
    }

else:
    update_experiment = {
        'fields': update_dict
    }


# Send the updates
req.put(added_exp__url, headers=auth_header, json=update_experiment)
