import ssl

from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import os
import argparse
ssl._create_default_https_context = ssl._create_unverified_context

PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here
import pandas as pd
# Set the BLAST parameters
blast_program = "blastp" # assuming you are working with amino acid sequences
blast_database = "nr" # replace with the actual name of your BLAST database
query_sequence = "mqrlvawdpa clplpppppa fksmevanfy yeadclaaay ggkaapaapp aarpgprppa gelgsigdhe raidfspyle plgapqapap atatdtfeaa ppapapapas sgqhhdflsd " \
                 "lfsddyggkn ckkpaeygyv slgrlgaakg alhpgcfapl hppppppppp aelkaepgfe padckrkeea gapgggagma agfpyalray lgyqavpsgs sgslstssss sppgtpspad akapptacya gaapapsqvk skakktvdkh sdeykirrer nniavrksrd kakmrnletq hkvleltaen erlqkkveql srelstlrnl fkqlpeplla ssghc"
query_coverage = 80 # replace with the desired query coverage percentage

# Run the BLAST search and retrieve the results
result_handle = NCBIWWW.qblast(program="blastp", database="nr", sequence=query_sequence, hitlist_size=100)
blast_records = NCBIXML.parse(result_handle)

# Iterate over the results and filter by query coverage
results= []
results_24 = []
results_199 = []
#TODO: change to 3 lists acording to hsp.query_start
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if (hsp.query_end - hsp.query_start + 1) / float(len(query_sequence)) >= query_coverage / 100.0:
                results.append((alignment.title, hsp.score))

# Sort the results by score
results = sorted(results, key=lambda x: x[1], reverse=True)

# Retrieve the top 5000 results
top_results = results[:5000]
print(top_results)
#TODO: put each list in diffrent file

with open("blast_results.txt ") as f:
    for title, score in top_results:
        print(title)
        f.write(title)
    f.close()

