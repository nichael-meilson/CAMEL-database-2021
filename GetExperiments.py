# Script to take data from API and export to Excel so I can get relevant information rather than hand searching through
# the database
import json
import requests
import pandas as pd

experiments = requests.get("https://cameldatabase.com/CAMEL/api/experiment?1=Saccharomyces%20cerevisiae").json()
eco = requests.get("https://cameldatabase.com/CAMEL/api/experiment?1=Escherichia%20coli").json()

df1 = pd.DataFrame.from_dict(experiments)
df2 = pd.DataFrame.from_dict(eco)

frames = [df1, df2]
df = pd.concat(frames)
df11 = pd.DataFrame(df['name'])
df22 = pd.DataFrame(df['fields'].values.tolist())
df33 = pd.DataFrame(df['id'])
df44 = pd.DataFrame(df['references'].values.tolist())
df55 = pd.DataFrame(df44[0].values.tolist())
df_merge = df11.reset_index().merge(df22.reset_index(), left_index=True, right_index=True, how='left')
df_merge1 = df_merge.reset_index().merge(df33.reset_index(), left_index=True, right_index=True, how='left')
df_merge2 = df_merge1.reset_index().merge(df55.reset_index(), left_index=True, right_index=True, how='left')
print(df_merge2)
df_merge2.to_excel("Experiment_information" + '.xlsx', index=False)