#! /usr/bin/python
from Bio import Phylo
import phyloK2 as K2
import sys
import glob
import csv

print("Initializing...")

tree_dir = sys.argv[1]
paths = glob.glob(tree_dir + "*.tre")
tree_names = [s[len(tree_dir):-4] for s in paths]

trees = [Phylo.read(x, 'newick') for x in paths]

k2 = K2.PhyloKernel(scaleFactor=0.3, gaussFactor=2)

for tree in trees:
	tree.ladderize()

if sys.argv[2] != "":
	k = K2.PhyloKernel(scaleFactor=0.3, gaussFactor=2, labelFactor=1, labelFilter=sys.argv[2])
	print("Computing kernels (labelled)...")

	with open(tree_dir + "kernels.csv", 'w') as f:
		w = csv.writer(f)
	
		w.writerow(tree_names)
		
		for i in xrange(0, len(trees)):		
			kernels = [k.kernel(trees[i], trees[j]) if j >= i else 0 for j in xrange(0, len(trees))]

			row = [str(x) for x in kernels]
	
			w.writerow(row)
		
			if i % 10 == 9:
				print i + 1
				f.flush()
			
print("\nComputing kernels (unlabelled)...")

with open(tree_dir + "kernels.u.csv", 'w') as f:
	w = csv.writer(f)

	w.writerow(tree_names)
		
	for i in xrange(0, len(trees)):	
	
		kernels = [k2.kernel(trees[i], trees[j]) if j >= i else 0 for j in xrange(0, len(trees))]

		row = [str(x) for x in kernels]
	
		w.writerow(row)
		
		if i % 10 == 9:
			print i + 1
			f.flush()