from Bio import Phylo
import re
from StringIO import StringIO

pat = re.compile('\[\&([^,\]]+),')

handle = open('../data/RPTree_output.nexus', 'rU')
text = handle.read()
handle.close()

matches = re.findall(pat, text)
while matches:
    # search and replace
    text = re.sub(pat, '[&\\1 ', text)
    matches = re.findall(pat, text)

# check for balanced parentheses


try:
    tree = Phylo.read(StringIO(text), 'nexus')
except:
    raise

#tree = Phylo.read('RPTree_output.nexus', 'nexus')


