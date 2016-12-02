"""
Automated Genbank queries to grab collection dates when available
Requires lxml
"""
import os
import sys
from seqUtils import convert_fasta
import re
import gb
from time import sleep
from csv import DictWriter

# regex for matching accession numbers
p = re.compile('[bdegjm]+_([A-Z]+[0-9]+)|(NC_[0-9]+)|([A-Z]{1,2}[0-9]{5,6})[^0-9]|([A-Z]{1,2}[0-9]{5,6})$')

p1 = re.compile('Sample timepoint ([0-9]+-\w+-[0-9]+)')
p2 = re.compile('Sample tissue (\w+)')
p3 = re.compile('Patient code ([0-9]+)')

fieldnames = ['header', 'accno', 'patid', 'coldate', 'tissue', 'moltype', 'isolate', 'isosource']

infile = open(sys.argv[1], 'rU')
fasta = convert_fasta(infile.readlines())
infile.close()

outfile = open(sys.argv[2], 'w')
writer = DictWriter(outfile, fieldnames=fieldnames)
writer.writeheader()

for i, (h, s) in enumerate(fasta):
    matches = p.findall(h)
    if len(matches) == 0:
        print 'unable to find accession number', h
        continue
    
    acc = max(matches[0], key=len)
    print acc
    
    entry = gb.GenbankEntry(acc)
    if entry.xml is None:
        # failed to retrieve record
        outfile.write('%s,NA\n' % h)
        continue
    
    if not hasattr(entry, 'source_dict'):
        outfile.write('%s,,,,\n' % h)
        continue

    date = entry.source_dict.get('collection_date', '')
    country = entry.source_dict.get('country', '')
    host = entry.source_dict.get('host', '')
    mol_type = entry.source_dict.get('mol_type', '')
    isolate = entry.source_dict.get('isolate', '')
    isolation_source = entry.source_dict.get('isolation_source', '')
    group = entry.source_dict.get('group', '')
    
    match = p1.findall(entry.comment)
    coldate = match[0] if match else ''
    match = p2.findall(entry.comment)
    tissue = match[0] if match else ''
    match = p3.findall(entry.comment)
    patid = match[0] if match else ''
    
    writer.writerow({
        'header': h,
        'accno': acc,
        'patid': patid if patid else group,
        'coldate': coldate if coldate else date,
        'tissue': tissue,
        'moltype': mol_type,
        'isolate': isolate,
        'isosource': isolation_source
    })
    
    sleep(0.1) # wait a second, to avoid spamming Genbank

outfile.close()

