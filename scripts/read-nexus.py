from Bio import Phylo
import re
from StringIO import StringIO
import sys

pat = re.compile('\[\&([^,\]]+),')

handle = open('../data/RPTree_output.nexus', 'rU')
text = handle.read()
handle.close()

matches = re.findall(pat, text)
while matches:
    # search and replace
    text = re.sub(pat, '[&\\1 ', text)
    matches = re.findall(pat, text)


tree = Phylo.read(StringIO(text), 'nexus')

