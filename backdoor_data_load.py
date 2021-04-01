## Input: 

from Addexperiment import _get_all_mutation_data, get_data_and_add_experiment
import pandas as pd
import numpy as np
import sqlalchemy

def rename_keys(dict_row):
    new_key = {
            '_1': 'Start POS',
            '_2': 'End POS',
            '_7': '∆AA',
            '_27': 'Category of Mutation',
            '_28': 'Interactions/Score',
            '_29': 'Total Interaction Score',
            '_30': 'Cello2go probabilities'
        }

    for item in list(dict_row.keys()):
        if item.startswith('_'):
            dict_row[new_key[item]] = dict_row.pop(item)
    
    return dict_row


def connect_to_db():
    config = {
    'host': 'localhost',
    'port': 3306,
    'user': 'camel',
    'password': 'abcdef',
    'database': 'CAMEL'
    }
    db_user = config.get('user')
    db_pwd = config.get('password')
    db_host = config.get('host')
    db_port = config.get('port')
    db_name = config.get('database')

    connection_str = f'mysql+pymysql://{db_user}:{db_pwd}@{db_host}:{db_port}/{db_name}'

    engine = sqlalchemy.create_engine(connection_str)
    return engine


def load_mutation_data(mutfile, exp_id):
    # Get data from xlsx files
    # run through mutation functions
    df = _get_all_mutation_data(mutfile)

    # take only columns needed
    field_mapper = {
            'CHROM_ID': (41, 'value_VARCHAR'),
            'Start POS': (42, 'value_INT'),
            'End POS': (43, 'value_INT'),
            'TYPE': (44, 'value_VARCHAR'),
            'REF': (45, 'value_TEXT'),
            'ALT': (46, 'value_TEXT'),
            'GEN': (47, 'value_TEXT'),
            '∆AA': (48, 'value_VARCHAR'),
            'POP': (49, 'value_VARCHAR'),
            'CLON': (50, 'value_VARCHAR'),
            'TIME': (51, 'value_VARCHAR'),
            'FREQ': (52, 'value_VARCHAR'),
            'COM': (53, 'value_BOOL'),
            'impact': (54, 'value_DOUBLE'),
            'score': (55, 'value_DOUBLE'),
            'ic': (56, 'value_DOUBLE'),
            'ddg': (57, 'value_DOUBLE'),
            'pdb_id': (58, 'value_VARCHAR'),
            'sequence': (59, 'value_VARCHAR'),
            'accession': (60, 'value_VARCHAR'),
            'modification_type': (61, 'value_VARCHAR'),
            'site_function': (62, 'value_VARCHAR'),
            'function_evidence': (63, 'value_VARCHAR'),
            'predicted_kinase': (64, 'value_VARCHAR'),
            'probability_loss': (65, 'value_DOUBLE'),
            'knockout_pvalue': (66, 'value_DOUBLE'),
            'tf': (67, 'value_VARCHAR'),
            'Category of Mutation': (68, 'value_TEXT'),
            'Interactions/Score': (69, 'value_DOUBLE'),
            'Total Interaction Score': (70, 'value_DOUBLE'),
            'Cello2go probabilities': (71, 'value_VARCHAR'),
            'Location': (72, 'value_VARCHAR')
        }

    field_df = df[field_mapper.keys()]

    # id: auto incremented?
    # experiment_id: do some querying to manually increment

    field_columns = ['experiment_id','mutation_id','field_id','value_VARCHAR','value_TEXT','value_INT','value_DOUBLE','value_BOOL','value_ATTACH']
    formatted_df = pd.DataFrame(columns=field_columns)

    # transform data to mutations_fields format
    import datetime
    mut_id = 1
    for rows in field_df.itertuples(index=False):
        print(datetime.datetime.now(), mut_id, "/", len(field_df))
        dict_row = rows._asdict()
        dict_row = rename_keys(dict_row)
        for key in dict_row.keys():
            tmp_list = [np.nan]*len(field_columns)
            index = field_columns.index(field_mapper[key][1])
            tmp_list[0] = exp_id
            tmp_list[1] = mut_id
            tmp_list[2] = field_mapper[key][0]
            tmp_list[index] = dict_row[key]
            tmp_series = pd.Series(tmp_list, index=formatted_df.columns)
            formatted_df = formatted_df.append(tmp_series, ignore_index=True)
        mut_id += 1

    engine = connect_to_db()

    formatted_df.to_sql('mutations_fields', con=engine, if_exists="replace")


def load_experiment(metafile, mutfile):
    answer = get_data_and_add_experiment(metafile)
    exp_id = answer['id']
    load_mutation_data(mutfile, exp_id)