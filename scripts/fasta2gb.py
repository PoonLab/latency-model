"""
Convert FASTA into a list of genbank accession numbers
"""
from seqUtils import convert_fasta
import sys
import re

pat = re.compile('[A-Z]{2}[0-9]{6}')

infile = open(sys.argv[1], 'rU')
outfile = open(sys.argv[2], 'w')


fasta = convert_fasta(infile)

for h, s in fasta:
    match = pat.findall(h)
    if match:
        outfile.write('%s\n' % (match[0], ))

outfile.close()
