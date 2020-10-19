# CAMEL-database-
Thesis scripts for 2019-2020 Bioinformatics Masters

Consists of the five Python scripts used throughout the thesis as well as a Jupyter Notebook where some analysis is done on the collected data.

## AddExperiment.py 
This is the wrap script that needs to run to work with CAMEL. Requires three other scripts to perform its functionality (CellularLocation.py, Mechismo_functionality.py, and MutFunc_functionality.py)

Consits of multiple functions to work is used to work with the CAMEL API and add/delete experiments and mutation data. The file consists of three main functions:

  remove_experiment(eid) - Removes an experiment from CAMEL, where eid is the experiment ID
  
  get_data_and_add_experiment(file, mutfile) - Takes in an experiment metadata file (file) and an optional mutation file (mutfile) and creates a new entry on CAMEL with the given metadata as well as uploading the mutation file and a new updated mutation file with analysis done on the mutations by CELLO2GO, mutfunc and mechismo.
  
  add_mutation_to_experiment(mutation_file) - Given a mutation file that is named appropriately (experimentID_filename) will upload the mutation file to the matching experiment on CAMEL and perform analysis on the file by CELLO2GO, mutfunc and mechismo.

## CellularLocation.py
Consists of two functions that take in a set of mutations in the CAMEL format and submit them one by one to CELLO2GO and return the respective subcellular locations of the targeted mutations. Can be run separately but the main design is to run in AddExperiment.py as part of adding mutation data to an experiment or creating an experiment with mutation data.

  locations(file) - Takes in an excel mutation file (file) in the CAMEL format and returns a list of the fasta sequences for each mutation.
  
  cello2go(genes) - Goes through the given list of fasta sequences (genes) that is outputted by the locations function and collects the subcellular location one by one running the input through CELLO2GO
  
## MutFunc_functionality.py
Consists of three functions that take in a set of mutations in the CAMEL format and submit them to mutfunc to get the ouput for these mutations. Can be run separately but the main design is to run in AddExperiment.py as part of adding mutation data to an experiment or creating an experiment with mutation data.

  runmutfunc(file) - Takes in a excel file of mutations in the CAMEL format (file) and returns the output of the mutfunc results
  
  extract_files(mut_func_file, mutation_data_frame) - takes in the ouput of the mutfunc results (mut_func_file) and the dataframe created from the file of mutations and extracts the results of mutfunc and places them in the dataframe.
  
  add_column_description() - Adds comments to all the new columns created in MutFunc_functionality.py, CellularLocation.py and Mechismo_functionality.py to describe what information is in these new columns.
  
## Mechismo_functionality.py
  Consists of one function that takes in a set of mutations in the CAMEL format and submits it to mechismo for the output for these mutations. Meant to be combined with CellularLocation.py and MutFunc_functionality.py in the AddExperiment.py scripts.
  
    run_mechismo(mutation_file) - takes in a excel mutation file in the CAMEL format (mutation_file) and outputs a dataframe that is the original mutation information alongside the mechismo output.
  
## AutomateMutationDAta.py
A script used to convert the mutation data from ALEdb to the desired CAMEL output in excel.
  
    Consists of two functions:
    
      The first mutTranslate(csvfile) - Takes in a mutation file taken from ALEdb as a csv file (csvfile) and converts it to an excel file in the CAMEL mutation format.
      
      append_df_to_excel() - Appends a dataframe into an existing excel file, it is used to add the converted mutations to the CAMEL mutation file with the existing template.
      
The final file of interest is the DataAnalysis.ipynb where all of the data analysis was done throughout this thesis.

  
