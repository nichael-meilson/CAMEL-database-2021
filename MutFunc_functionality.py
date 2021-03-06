# Goal is to take a file of mutation data, retrieve all the unique SNPs from it and run them through
# MutFunc(http://www.mutfunc.com/ to provide extra information to the data providers about their mutations
# This script is ideally attached to the button the experiment page that will return fresh data to the user on demand
import pandas as pd
import mechanize
import time
import urllib.request
import zipfile
import os
from Bio.Data.IUPACData import protein_letters_3to1
# import win32com.client


# requires an excel file of mutations that should be obtained from the Camel database experiment page. Since MutFunc
# only works on SNPs the first thing we have to do is remove all non SNPs and only keep unique values

def extract_snps(df):
    # Input: mutfile df
    # Output: filted mutfile df for SNPs
    SNP_df = df[(df['TYPE']=="SNP") 
                & (df["GEN"]!="NA")
                & (df["GEN"].notnull())]
    SNP_df.drop_duplicates(inplace=True)
    return(SNP_df)

def convert_snp_mutfunc(df):
    # Convert SNPs to a suitable input format for Mutfunc
    # Take only data from the SNPs that we want for checking in MutFunc
    # Split numbers and test in cell (ie. 'Pro375Ser' --> ['Pro','375','Ser'])
    df_split = df["∆AA"].str.findall(r'[A-Za-z]+|\d+')

    # Iterate through series to add empty values to get matching array lengths 
    # ie. ['Lys','8'] --> ['Lys','8','']
    # ie. ['Pro'] --> ['Pro','','']
    # Results in a DF with "reference","position","mutant" as columns
    indexed_list = list()
    for entry in df_split:
        if len(entry) == 2:
            entry.append("")
        elif len(entry) == 1:
            entry.append("")
            entry.append("")
        indexed_list.append(entry)

    aa_mutant_df = pd.DataFrame(indexed_list, columns=["reference","position","mutant"])

    # Replace 3 letter IUPAC protein codes with 1 letter protein codes
    aa_mutant_df.replace({"reference": protein_letters_3to1, "mutant": protein_letters_3to1}, inplace=True)

    # Append chromosomes to this DF to determine which genes don't have any specific mutation data
    aa_mutant_df["GEN"] = df["GEN"]
    final_df = aa_mutant_df[(aa_mutant_df["mutant"].str.len() == 1) 
                            & (aa_mutant_df["GEN"].notnull())
                            & (aa_mutant_df["GEN"] != "")]
    SNP_list = pd.DataFrame()
    # Merge all parts together so they can be added as one entry per line
    SNP_list["SNPs"] = final_df["GEN"] + " " + final_df["reference"] + final_df["position"] + final_df["mutant"]
    
    # I do not think we need to keep duplicated entries here as we add back in on start pos in the second function
    SNP_list = SNP_list.drop_duplicates()
    return(SNP_list)


def runmutfunc(file):
    df = pd.read_excel(file, header=4, keep_default_na=False)
    # Here we filter to only keep SNP mutations
    snps_df = extract_snps(df)
    SNP_list = convert_snp_mutfunc(snps_df)

    # time now to learn how to access the website through python, looking at the mechanize package
    br = mechanize.Browser()
    br.open("http://www.mutfunc.com/submit#")
    br.select_form(nr=0)
    species = br.form.find_control("tax_id")
    for item in species.items:
        if item.name == "2":
            item.selected = True
    mutations = ""
    for i in range(len(SNP_list)):
        mutations = mutations + SNP_list["SNPs"].iloc[i] + "\n"
    br.form['muts_text'] = mutations

    # If the number of mutations is too large this will not work. A way of handling this would be to split your
    # mutations into parts and then later on go through each part of the results you save. Needs to be extremely large
    # though as 1500+ mutations works.

    br.submit()
    base_url = br.geturl()
    # SNP_list.to_csv('test.csv',index=False, header=False)
    # cannot figure out how to have the url updated from the wait page to the results page, even though it is redirected
    # so we will just set an arbitrarily long wait time where we then assume that the loading as finished and we move to
    # the export page which lets us download the files
    # Here we pause for 3 minutes (can be changed) before switching our page to the results
    time.sleep(180)
    new_url = str(base_url).replace("wait", "export")
    urllib.request.urlretrieve(new_url, file + ".gz")
    mut_func_file = file + ".gz"
    return extract_files(mut_func_file, df)
    # Here we then want to save this file so that it can be added to the field when experiments are added or mutation
    # data is attached
    # extract_files("C:/Users/samue/Desktop/Thesis/35_42C.csv.xlsx.gz", mutation_update_df)
    # this part on is for testing and doesnt need to necessarily be part of the add experiment though maybe it could
    # updated_data_frame.to_excel("Mutations with Mutfunc results.xlsx", index=False)
    # zip = zipfile.ZipFile("C:/Users/samue/Desktop/Thesis/35_42C.csv.xlsx.gz", 'a')
    # zip.write("Mutations with Mutfunc results.xlsx")
    # zip.close()

def convert_aa_column_3to1letter(df, column='∆AA'):
    # Converts the given column of a dataframe from having 3 letter IUPAC protein codes to 1 letter codes
    # Returns the same input df but with the converted '∆AA' column
    # Ex: Phe455Val --> F445V

    # Mutfunc switched to 1 letter codes
    # Split numbers and test in cell (ie. 'Pro375Ser' --> ['Pro','375','Ser'])
    df_split = df[column].str.findall(r'[A-Za-z]+|\d+')

    # Iterate through series to add empty values to get matching array lengths 
    # ie. ['Lys','8'] --> ['Lys','8','']
    # ie. ['Pro'] --> ['Pro','','']
    # Results in a DF with "reference","position","mutant" as columns
    indexed_list = list()
    for entry in df_split:
        if len(entry) == 2:
            entry.append("")
        elif len(entry) == 1:
            entry.append("")
            entry.append("")
        indexed_list.append(entry)

    aa_mutant_df = pd.DataFrame(indexed_list).loc[:,:2]
    aa_mutant_df.columns = ["reference","position","mutant"]

    # Replace 3 letter IUPAC protein codes with 1 letter protein codes
    # aa_mutant_df.replace({"reference": protein_letters_3to1, "mutant": protein_letters_3to1}, inplace=True)
    # aa_mutant_df
    for i in range(0,len(aa_mutant_df)):
        try:
            aa_mutant_df.iloc[i]["reference"] = protein_letters_3to1[aa_mutant_df.iloc[i]["reference"]]
            aa_mutant_df.iloc[i]["mutant"] = protein_letters_3to1[aa_mutant_df.iloc[i]["mutant"]]
        except:
            pass

    df[column] = aa_mutant_df["reference"]+aa_mutant_df["position"]+aa_mutant_df["mutant"]
    df["single_letter_reference"] = aa_mutant_df["reference"]
    df["single_letter_position"] = aa_mutant_df["position"]
    df["single_letter_mutant"] = aa_mutant_df["mutant"]

    return(df)


def extract_files(mut_func_file, mutation_data_frame):
    # now we have the zipped file and we want to unzip it and reach each file one by one and collect the important parts
    # we know the order of the unzipped files which is psites, start_stop, interfaces, other_ptms, linear_motifs,
    # conservation, stability, tfbs
    # add new columns to our data frame that will be our final results
    # need to also add information about what each column represents
    df = mutation_data_frame
    # Remove intergenic mutations, but keep everything else
    df = df[df["GEN"] != 'NA']
    columns_to_add = ["refaa", "altaa", "impact", "score", "ic", "ddg", "pdb_id", "sequence", "accession",
                      "modification_type", "site_function", "function_evidence", "predicted_kinase", "probability_loss",
                      "knockout_pvalue", "tf", "Category of Mutation"]
    df = df.reindex(columns=df.columns.tolist() + columns_to_add)
    df = convert_aa_column_3to1letter(df)
    df.fillna('', inplace=True)

    zip_file_object = zipfile.ZipFile(mut_func_file, 'r')
    # Each file has a different header length so we will do each individually as well as different requirements of what
    # data to retrieve
    # first we do psites
    p_sites = zip_file_object.open(zip_file_object.namelist()[0])
    p_sites_mutations = pd.read_csv(p_sites, skiprows=24, header=0, delimiter='\t')
    # check to see if the file is empty outside of the header
    if p_sites_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(p_sites_mutations)):
                if df.loc[i, "GEN"] == p_sites_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == p_sites_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == p_sites_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = p_sites_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = p_sites_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = p_sites_mutations.loc[j, "impact"]
                    df.loc[i, "site_function"] = p_sites_mutations.loc[j, "site_function"]
                    df.loc[i, "function_evidence"] = p_sites_mutations.loc[j, "function_evidence"]
                    df.loc[i, "predicted_kinase"] = p_sites_mutations.loc[j, "predicted_kinase"]
                    df.loc[i, "probability_loss"] = p_sites_mutations.loc[j, "probability_loss"]
                    if df.loc[i, "Category of Mutation"] != "" and "Psites" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Psites"
                    else:
                        df.loc[i, "Category of Mutation"] = "Psites"
    p_sites.close()

    # now we do start_stop
    start_stop = zip_file_object.open(zip_file_object.namelist()[1])
    # Need to try/except block this file because it is the only one that does not include a header for Mutfunc results
    # and can cause some issues because of this since it has been empty most of the time
    try:
        start_stop_mutations = pd.read_csv(start_stop, skiprows=11, header=None, delimiter='\t', index_col=False)
        if start_stop_mutations.empty:
            pass
        else:
            for i in range(0, len(df)):
                for j in range(0, len(start_stop_mutations)):
                    if df.loc[i, "GEN"] == start_stop_mutations.loc[j, "gene"]\
                        and df.loc[i, "single_letter_reference"] == start_stop_mutations.loc[j, "refaa"]\
                        and df.loc[i, "single_letter_mutant"] == start_stop_mutations.loc[j, "altaa"]:
                        df.loc[i, "refaa"] = start_stop_mutations.loc[j, "refaa"]
                        df.loc[i, "altaa"] = start_stop_mutations.loc[j, "altaa"]
                        df.loc[i, "impact"] = start_stop_mutations.loc[j, "impact"]
                        if df.loc[i, "Category of Mutation"] != "" and "Start_Stop" not in df.loc[i,'Category of Mutation']:
                            df.loc[i, "Category of Mutation"] += ", Start_Stop"
                        else:
                            df.loc[i, "Category of Mutation"] = "Start_Stop"
    except:
        pass
    start_stop.close()

    interfaces = zip_file_object.open(zip_file_object.namelist()[2])
    interfaces_mutations = pd.read_csv(interfaces, skiprows=20, header=0, delimiter="\t", index_col=False)
    if interfaces_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(interfaces_mutations)):
                if df.loc[i, "GEN"] == interfaces_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == interfaces_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == interfaces_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = interfaces_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = interfaces_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = interfaces_mutations.loc[j, "impact"]
                    df.loc[i, "pdb_id"] = interfaces_mutations.loc[j, "pdb_id"]
                    df.loc[i, "ddg"] = interfaces_mutations.loc[j, "ddg"]
                    if df.loc[i, "Category of Mutation"] != "" and "Interfaces" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Interfaces"
                    else:
                        df.loc[i, "Category of Mutation"] = "Interfaces"
    interfaces.close()

    other_ptms = zip_file_object.open(zip_file_object.namelist()[3])
    other_ptms_mutations = pd.read_csv(other_ptms, skiprows=16, header=0, delimiter="\t", index_col=False)
    if other_ptms_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(other_ptms_mutations)):
                if df.loc[i, "GEN"] == other_ptms_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == other_ptms_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == other_ptms_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = other_ptms_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = other_ptms_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = other_ptms_mutations.loc[j, "impact"]
                    df.loc[i, "modification_type"] = other_ptms_mutations.loc[j, "modification_type"]
                    if df.loc[i, "Category of Mutation"] != "" and "Other_PTM" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Other_PTM"
                    else:
                        df.loc[i, "Category of Mutation"] = "Other_PTM"
    other_ptms.close()

    linear_motifs = zip_file_object.open(zip_file_object.namelist()[4])
    linear_motifs_mutations = pd.read_csv(linear_motifs, skiprows=21, header=0, delimiter="\t", index_col=False)
    if linear_motifs_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(linear_motifs_mutations)):
                if df.loc[i, "GEN"] == linear_motifs_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == linear_motifs_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == linear_motifs_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = linear_motifs_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = linear_motifs_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = linear_motifs_mutations.loc[j, "impact"]
                    df.loc[i, "sequence"] = linear_motifs_mutations.loc[j, "sequence"]
                    df.loc[i, "accession"] = linear_motifs_mutations.loc[j, "accession"]
                    if df.loc[i, "Category of Mutation"] != "" and "Linear_Motif" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Linear_Motif"
                    else:
                        df.loc[i, "Category of Mutation"] = "Linear_Motif"
    linear_motifs.close()

    conservation = zip_file_object.open(zip_file_object.namelist()[5])
    conservation_mutations = pd.read_csv(conservation, skiprows=16, header=0, delimiter="\t", index_col=False)
    if conservation_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(conservation_mutations)):
                if df.loc[i, "GEN"] == conservation_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == conservation_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == conservation_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = conservation_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = conservation_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = conservation_mutations.loc[j, "impact"]
                    df.loc[i, "score"] = conservation_mutations.loc[j, "score"]
                    df.loc[i, "ic"] = conservation_mutations.loc[j, "ic"]
                    if df.loc[i, "Category of Mutation"] != "" and "Conservation" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Conservation"
                    else:
                        df.loc[i, "Category of Mutation"] = "Conservation"
    conservation.close()

    stability = zip_file_object.open(zip_file_object.namelist()[6])
    stability_mutations = pd.read_csv(stability, skiprows=17, header=0, delimiter="\t", index_col=False)
    if stability_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(stability_mutations)):
                if df.loc[i, "GEN"] == stability_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == stability_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == stability_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = stability_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = stability_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = stability_mutations.loc[j, "impact"]
                    df.loc[i, "pdb_id"] = stability_mutations.loc[j, "pdb_id"]
                    df.loc[i, "ddg"] = stability_mutations.loc[j, "ddg"]
                    if df.loc[i, "Category of Mutation"] != "" and "Stability" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", Stability"
                    else:
                        df.loc[i, "Category of Mutation"] = "Stability"
    stability.close()

    tfbs = zip_file_object.open(zip_file_object.namelist()[7])
    tfbs_mutations = pd.read_csv(tfbs, skiprows=28, header=0, delimiter="\t", index_col=False)
    if tfbs_mutations.empty:
        pass
    else:
        for i in range(0, len(df)):
            for j in range(0, len(tfbs_mutations)):
                if df.loc[i, "GEN"] == tfbs_mutations.loc[j, "gene"]\
                    and df.loc[i, "single_letter_reference"] == tfbs_mutations.loc[j, "refaa"]\
                    and df.loc[i, "single_letter_mutant"] == tfbs_mutations.loc[j, "altaa"]:
                    df.loc[i, "refaa"] = tfbs_mutations.loc[j, "refaa"]
                    df.loc[i, "altaa"] = tfbs_mutations.loc[j, "altaa"]
                    df.loc[i, "impact"] = tfbs_mutations.loc[j, "impact"]
                    df.loc[i, "tf"] = tfbs_mutations.loc[j, "tf"]
                    df.loc[i, "knockout_pvalue"] = tfbs_mutations.loc[j, "knockout_pvalue"]
                    if df.loc[i, "Category of Mutation"] != "" and "TFBS" not in df.loc[i,'Category of Mutation']:
                        df.loc[i, "Category of Mutation"] += ", TFBS"
                    else:
                        df.loc[i, "Category of Mutation"] = "TFBS"
    tfbs.close()

    # Drop the columns in the df used for the matching/indexing
    df.drop(columns=["single_letter_reference", "single_letter_mutant", "single_letter_position"], inplace=True)
    return df


def add_column_description():

    xl = win32com.client.Dispatch("Excel.Application")
    xl.Visible = 1
    current_file_path = os.getcwd() + "\Mutation_results.xlsx"
    wb = xl.Workbooks.Open(current_file_path)
    sheet = wb.ActiveSheet
    # add comments
    sheets = ["O1", "P1", "Q1", "R1", "S1", "T1", "U1", "V1", "W1", "X1", "Y1", "Z1", "AA1",
              "AB1", "AC1", "AD1", "AE1", "AF1", "AG1", "AH1", "AI1"]
    comments = ["Reference amino acid", "Mutated amino acid",
                "Is the mutation predicted to impact function? '1' if yes, '0' if no",
                "Sift score, any mutation with a score below 0.05 is considered deleterious ",
                "Information content at this position of the alignment (a high value indicates strong conservation, where the maximum value is 4.32)",
                "Predicted change in free energy of unfolding, where a value above 0 indicates a destabilising mutation",
                "Pdb identifier or homology model identifier of the structure containing the mutation",
                "Sequence of the linear motif", "ELM accession for the linear motif",
                "Type of post-translational modifications", "Function of this phosphorylation site, if any",
                "Evidence of site function, if any", "Kinase predicted to lose phosphorylation at this site",
                "Probability of kinase losing phosphorylation at this site",
                "P-value of over or under-expression for the downstream gene when the transcription factor is knocked out",
                "Transcription factor predicted to bind this binding site", "Category of mutfunc mutation",
                "protein-protein interactions: gene name of interactor, protein-chemical interactions: '[CHEM:type:id]'"
                "DNA/RNA interactions: '[DNA/RNA]', Mechismo score for the interaction - The higher the Mechismo Score,"
                " the more likely a particular mutation or modification is to affect interactions with other molecules."
                , "Total combined Mechismo interaction score for all molecules",
                "Likelihood for each location that mutation is there", "Most probable location"]
    for column, comment in zip(sheets, comments):
        sheet.Range(column).AddComment()
        sheet.Range(column).Comment.Visible = False
        sheet.Range(column).Comment.Text(comment)
    updated_file_path = os.getcwd() + "\Mutation_results_complete.xlsx"
    wb.SaveAs(updated_file_path)
    wb.Close()
    xl.Quit()


# runmutfunc("C:/Users/samue/Desktop/Thesis/ALEDB_conversion/Experiment_Data/42C.csv.xlsx")
# extract_files("C:/Users/samue/Desktop/Thesis/35_42C.csv.xlsx.gz",
#               "C:/Users/samue/Desktop/Thesis/ALEDB_conversion/Experiment_Data/42C.csv.xlsx")