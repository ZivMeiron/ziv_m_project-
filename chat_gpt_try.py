import ssl
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import numpy as np


ssl._create_default_https_context = ssl._create_unverified_context
PATH = os.getcwd()
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here

# Set the BLAST parameters
blast_program = "blastp" # assuming you are working with amino acid sequences
blast_database = "nr" # replace with the actual name of your BLAST database
query_sequence = "mqrlvawdpa clplpppppa fksmevanfy yeadclaaay ggkaapaapp aarpgprppa gelgsigdhe raidfspyle plgapqapap atatdtfeaa ppapapapas sgqhhdflsd " \
                 "lfsddyggkn ckkpaeygyv slgrlgaakg alhpgcfapl hppppppppp aelkaepgfe padckrkeea gapgggagma agfpyalray lgyqavpsgs sgslstssss sppgtpspad akapptacya gaapapsqvk skakktvdkh sdeykirrer nniavrksrd kakmrnletq hkvleltaen erlqkkveql srelstlrnl fkqlpeplla ssghc"
query_coverage = 80 # replace with the desired query coverage percentage

START_QUERY_ONE = 24
START_QUERY_TWO = 199

def get_blast_records(records_num):
    # Run the BLAST search and retrieve the results
    print("Trying to fetch results from blast.....")
    result_handle = NCBIWWW.qblast(program="blastp", database="nr", sequence=query_sequence, hitlist_size=records_num)
    blast_records = NCBIXML.parse(result_handle)
    print("Finshed  fetch results from blast.....")
    print(f"Fetchted {len(blast_records)} results")
    

    return blast_records

def parse_results(blast_records, query_coverage: int = query_coverage):
    """
    Parses the BLAST results and returns the top results for the given query coverage
    :param blast_records: BLAST results
    :param query_coverage: Query coverage
    :return: List of top results (query start, subject title, score)
    
    """
    # Iterate over the results and filter by query coverage
    results= []

    #TODO: change to 3 lists acording to hsp.query_start
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if (hsp.query_end - hsp.query_start + 1) / float(len(query_sequence)) >= query_coverage / 100.0:
                     # if hsp.expect < np.exp(-10)^:    
                    results.append((hsp.query_start,alignment.title, hsp.score))

    # Sort the results by score
    # results = sorted(results, key=lambda x: x[2], reverse=True)
    # top_results = results[:top_res_num]
    print(f"Found {len(results)} results")
    return results

def splits_results(results):
    results_query_one = []
    results_query_two = []
    results_query_three = []

    for query_start, title, score in results:
        if query_start < START_QUERY_ONE:
            results_query_one.append((title, query_start, score))
        elif query_start < START_QUERY_TWO:
            results_query_two.append((title, query_start, score))
        else:
            results_query_three.append((title, query_start, score))

    return results_query_one, results_query_two , results_query_three


def write_results(results, file_name):
    with open(file_name, 'w') as f:
        for title, query_start, score in results:
            f.write(f"{title} {query_start} {score}\n")



if __name__ == '__main__':
    records = get_blast_records(100)
    results = parse_results(records, query_coverage)
    results_query_one, results_query_two , results_query_three = splits_results(results)
    write_results(results_query_one, 'results_query_one.txt')
    write_results(results_query_two, 'results_query_two.txt')
    write_results(results_query_three, 'results_query_three.txt')