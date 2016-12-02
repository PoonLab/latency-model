"""
Partition data sets by host, annotate with collection dates
"""

from seqUtils import convert_fasta
from csv import DictReader

handle = open('../data/25712966/env-cellular.fa', 'rU')
fasta = convert_fasta(handle)
handle.close()

handle = open('../data/25712966/env-cellular.gb.csv', 'rU')
rows = DictReader(handle)
for row in rows:
    print row['accno']
    break

    