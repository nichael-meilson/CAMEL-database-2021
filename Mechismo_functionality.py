# Making Mechismo work automatically with CAMEL
# Mechismo uses ids for each session so we can try using mechanize like in MutFunc

import mechanize
import pandas as pd
import time
import os
import mutfunc_functionality

def _format_SNP_list(SNP_list):
    SNP_list[['gene', 'aa', 'pos']] = SNP_list['SNPs'].str.split(' ',expand=True)
    mutfunc_functionality.convert_aa_column_3to1letter(SNP_list, column='aa')
    series = SNP_list['gene'] + "/" + SNP_list['single_letter_reference'] + SNP_list['single_letter_position'] + SNP_list['single_letter_mutant']
    result = pd.DataFrame(series, columns=["SNPs"])
    result.fillna("", inplace=True)
    result['Start POS'] = SNP_list['pos']
    return result

# Mechismo works with SNPs or with PTMs so we take only the SNPs in our mutations
def run_mechismo(mutation_file):
    # Making Mechismo work automatically with CAMEL
    # Mechismo uses ids for each session so we can try using mechanize like in MutFunc

    import mechanize
    import pandas as pd
    import time
    import os

    # Mechismo works with SNPs or with PTMs so we take only the SNPs in our mutations

    df = pd.read_excel(mutation_file, header=4, keep_default_na=False)
    # Here we filter to only keep unique SNP mutations
    df = df[df['TYPE'] == "SNP"].drop_duplicates()
    # Also have to remove NAs in gene
    df = df[df['GEN'] != "NA"]
    # Take only data from the SNPs that we want for checking in MutFunc
    # Merge all parts together so they can be added as one entry per line
    SNP_list = pd.DataFrame()
    SNP_list["SNPs"] = df['GEN'] + " " + df.iloc[:, 7] + " " + df["Start POS"].astype(str)
    # Have to remove noncoding and pseudogene
    SNP_list= SNP_list[~SNP_list["SNPs"].str.contains('pseudogene')]
    SNP_list = SNP_list[~SNP_list["SNPs"].str.contains('noncoding')]
    SNP_list = SNP_list.drop_duplicates()
    
    SNP_list = _format_SNP_list(SNP_list)

    # Now that we have our list of mutations that match Mechismos guidelines we can start running the website
    br = mechanize.Browser()
    br.open("http://mechismo.russelllab.org/")
    # for form in br.forms():
    #     print(form)
    br.select_form(nr=0)
    mutations = ""
    for i in range(len(SNP_list)):
        mutations = mutations + SNP_list["SNPs"].iloc[i] + "\n"

    br.form["search"] = mutations
    br.form["search_name"] = "Experiment Test"
    br.form["taxon"] = ["83333"]
    br.form["stringency"] = ["low"]
    br.submit()

    # No direct way to check its done so just pause x amount of time
    time.sleep(60)
    file_url = br.geturl()
    url = ""
    for i in file_url:
        url += i
    file_id = url.split("/")[-1]

    new_url = "http://mechismo.russelllab.org/static/data/jobs/" + str(file_id) + "/" + str(file_id) + ".site_table.tsv.gz"
    br.retrieve(new_url, "Experiment.tsv")

    results = pd.read_csv("Experiment.tsv", sep="\t", header=0)
    final_results = results[["name_a1", "name_b1", "mechProt", "mechChem", "mechDNA/RNA", "mech", "pos_a1"]]
    # Want to remove NaN from interactors to just get relevant information
    interactors = final_results['name_b1'].notna()
    updated_results = final_results[interactors]
    updated_results = updated_results.drop_duplicates()
    updated_results = updated_results.reset_index().drop(columns='index')

    # Remove [PROT] results
    updated_results = updated_results[updated_results['name_b1'] != "[PROT]"]
    camel_results = pd.DataFrame(columns=['Position', "Interactions/Score", "Total Interaction Score"])

    index = 0
    input_tracker = ''
    input_value = ""
    total_score = ""
    for index, row in updated_results.iterrows():
        if not input_tracker == '':
            if row["pos_a1"] == input_tracker:
                if "CHEM" in row["name_b1"]:
                    input_value += ", " + row["name_b1"] + ", " + str(row["mechChem"])
                elif "DNA" in row["name_b1"]:
                    input_value += ", " + row["name_b1"] + ", " + str(row["mechDNA/RNA"])
                else:
                    input_value += ", " + row["name_b1"] + ", " + str(row["mechProt"])
                continue
            else:
                camel_results.loc[index] = [input_tracker, input_value, total_score]
                index += 1
        input_tracker = row["pos_a1"]
        if "CHEM" in row["name_b1"]:
            input_value = row["name_b1"] + ", " + str(row["mechChem"])
        elif "DNA" in row["name_b1"]:
            input_value = row["name_b1"] + ", " + str(row["mechDNA/RNA"])
        else:
            input_value = row["name_b1"] + ", " + str(row["mechProt"])
        # This score could be a float instead
        total_score = str(row["mech"])
    # Have to add last entry
    camel_results.loc[index] = [input_tracker, input_value, total_score]
    # Have to manipulate dataframe to match mutation file so they can be merged by Start Pos
    camel_results = camel_results.rename(columns={'Position': "Start POS"})
    camel_results["Start POS"] = pd.to_numeric(camel_results["Start POS"], errors='ignore')
    # print(camel_results)
    return camel_results



