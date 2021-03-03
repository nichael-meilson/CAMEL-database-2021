# Will add a second source of data to researchers on CAMEL by providing information on the sub-cellular location of the
# genes that are being mutated, as well as functional gene ontology annotation. CELLO2GO
# Location of files needs to be updated throughout the script based on current directory used

import pandas as pd
import datetime
import re
import time
from selenium import webdriver
# from selenium.common.exceptions import ElementClickInterceptedException, NoSuchElementException
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait as wait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import pyperclip
import os
from zipfile import ZipFile
import pyautogui
import shutil
import concurrent.futures
import multiprocessing
import logging
import sys

# it is possible we can work with strains that are not the most popular since we just want fasta
def locations(file):
    df = pd.read_excel(file, header=4, keep_default_na=False)
    # The first thing to do is to remove all mutations that do not provide the gene that is affected as well as drop
    # duplicates
    genes = df['GEN'] != "NA"
    df = df[genes]
    df = df.drop_duplicates()
    df = df.reset_index()
    df = df.drop(columns="index")
    # only works for four strains downloaded
    strain_id = df['CHROM_ID'].unique()[0]
    if strain_id == "NC_000913":
        reference_genome = open("reference_genomes/EcoliNC000913.txt", 'r', encoding='UTF-8')
    elif strain_id == "NC_007779":
        reference_genome = open("reference_genomes/EscherichiacoliNC007779.txt", 'r', encoding='UTF-8')
    elif strain_id == "REL606":
        reference_genome = open("reference_genomes/EscherichiacoliBstr.REL606.txt", 'r', encoding='UTF-8')
    elif strain_id == "CP009273":
        reference_genome = open("reference_genomes/EColiCP009273.txt", 'r', encoding='UTF-8')
    # If none of these strains are in our dataset we end the function
    else:
        print("This strain of E. coli is not available right now")
        return "False"
    # creates dictionary for our reference strain with key gene name and value fasta sequence
    gene_dict = {}
    gene_name = None
    seq = ''
    take_seq = False
    for line in reference_genome:
        if line.startswith(">") and strain_id in line:
            if gene_name and not seq == '':
                gene_dict[gene_name.replace("-", "")] = seq
                seq = ''
            for item in line.split(" "):
                if 'gene=' in item:
                    gene_name = item.split('=')[1].strip("[]")
                    take_seq = True
        elif line.startswith(">") and strain_id not in line:
            take_seq = False
        else:
            if take_seq:
                seq += line.strip('\r\n')
    if gene_name and not seq == "":
        gene_dict[gene_name] = seq
    # collect all fasta sequences for each gene in our mutation list
    final_list = ""
    for rownumber, mutation in df.iterrows():
        for gene in re.split(', |;', mutation[6]):
            final_list += "> %s\n" % (mutation[1])
            try:
                final_list += str(gene_dict[gene]) + "\n"
            except KeyError:
                gene_synonyms = open("reference_genomes/Ecoli_gene_synonyms.tab", 'r', encoding='UTF-8')
                for line in gene_synonyms:
                    if gene in line:
                        names = line.split("\t")[1]
                        for name in names.split(" "):
                            try:
                                final_list += str(gene_dict[name]) + "\n"
                                break
                            except KeyError:
                                continue
                        break
                gene_synonyms.close()
    return final_list


# Now that we have the list of genes and their annotations we submit them to the cell2go website and get the resulting
# file, mechanize doesn't work so going to try selenium
# Depending on the browser you use you need to add an executable such as geckodriver (firefox) or chrome driver(chrome)
# to PATH, have to download file from https://github.com/mozilla/geckodriver/releases/tag/v0.26.0 (firefox)
# or https://sites.google.com/a/chromium.org/chromedriver/home (chrome) use
# this executable file in the path, Firefox wouldn't automatically save the file so switched to Chrome
# Throughout the function we use time.sleep to have the function wait for pages to load

def cello2go(genes):
    # Writes progress to a log file to see how far along the scraping is
    logging.basicConfig(filename='cello.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

    old_position = ""
    cell2go_columns = ["Start POS", "Cello2go probabilities", "Location"]
    location_results = pd.DataFrame(columns=cell2go_columns)
    updated = genes.split(">")[1:]
    total = len(updated)
    complete = 0
    logging.info('Starting scrape')
    for sequence in updated:
        try:
            sequence = ">" + sequence
            # Check to see if we did not find a gene name for our mutation and therefore we do not have a sequence
            if "G" not in sequence:
                continue
            position = (sequence.split("\n"))[0].strip(" ")
            if position == old_position:
                continue
            old_position = position
            # first we open up our webpage
            # This path needs to be where the chromedriver executable is stored
            path = "reference_genomes/chromedriver"
            options = webdriver.ChromeOptions()
            prefs = {
                "download.default_directory": r"reference_genomes\cello_files",
                "download.directory_upgrade": "true",
                "download.prompt_for_download": "false",
                "disable-popup-blocking": "true"
            }
            options.add_experimental_option("prefs", prefs)
            browser = webdriver.Chrome(executable_path=path, service_log_path='nul', options=options)
            browser.get("http://cello.life.nctu.edu.tw/cello2go/")
            # Then we clear the content and paste in our string of headers and FASTA sequences before running the search
            paste_sequence = browser.find_element_by_name("sequence")
            paste_sequence.clear()
            # Testing why one did not work
            # print(sequence)
            pyperclip.copy(sequence)
            print(sequence)
            paste_sequence.send_keys(Keys.COMMAND + "v")
            browser.find_element_by_xpath("/html/body/center/div[5]/form[1]/table/thead/tr/td[2]/button/span[2]").click()
            browser.find_element_by_xpath("/html/body/div[3]/ul/li[1]/label").click()
            submit_button = browser.find_element_by_id("do-blast")
            submit_button.click()
            time.sleep(9)
            location_values = []
            int_values = []
            # Scrape page
            ele = browser.find_elements_by_xpath("//table[@id='Bacteria-gramn']")
            for e in ele:
                for td in e.find_elements_by_xpath(".//td"):
                    location_values.append(td.text)
            # Update values from string to int to get max value
            for i in location_values[1::2]:
                int_values.append(float(i))
            # sometimes it doesn't pull the values so we can try this again if it fails then pull values a second time
            # just in case, testing this
            if not int_values:
                ele = browser.find_elements_by_xpath("//table[@id='Bacteria-gramn']")
                for e in ele:
                    for td in e.find_elements_by_xpath(".//td"):
                        location_values.append(td.text)
                # Update values from string to int to get max value
                for i in location_values[1::2]:
                    int_values.append(float(i))
            if max(int_values) == int_values[0]:
                mutation_location = "Extracellular"
            elif max(int_values) == int_values[1]:
                mutation_location = "Outermembrane"
            elif max(int_values) == int_values[2]:
                mutation_location = "Periplasmic"
            elif max(int_values) == int_values[3]:
                mutation_location = "Innermembrane"
            else:
                mutation_location = "Cytoplasmic"
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+": "+str(len(location_results)+1)+"/"+str(len(updated)))
            location_results.loc[len(location_results)] = [float(position.split(" ")[1]), ",".join(location_values), str(mutation_location)]
            browser.close()
            browser.quit()
        # just in case the sequence doesnt work this will then move on instead of stopping the program
        except:
            continue
        complete += 1
        logging.info('************************************************************************************************')
        logging.info(f'Complete Scrapes: {complete}/{total}')
    return location_results


def get_batch_list(genes, N=4):
    updated = genes.split(">")[1:]

    # For test purposes, use the first M genes
    # Comment this out later
    # M = 4
    # updated = updated[0:M]
    
    # We need the N input of batch_division to be the size of the batch, not the number of sequence batches
    total_batches = round(len(updated)/N)

    batch_list = list(batch_division(updated, total_batches))

    final_batch_list = []
    for sequences in batch_list:
        final_batch_list.append(">".join(sequences))

    final_batch_list = [">"+string for string in final_batch_list]

    return final_batch_list


def batch_division(l, n): 
    # This returns a generator - don't think too much about it
    for i in range(0, len(l), n):  
        yield l[i:i + n] 


def run_cello2go(file, N=4):
    # Runs multiple browsers for webscraping - cuts down data processing time significantly
    genes = locations("test/test_mutation_file.xlsx")
    batch_list = get_batch_list(genes, N)
    result_list = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(cello2go, batch_list)

        for result in results:
            result_list.append(result)

    df = pd.concat(result_list)
    df = df.reset_index().drop(columns='index')
    return df

