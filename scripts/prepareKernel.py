#! /usr/bin/python
from Bio import Phylo
import phyloK2 as K2
import sys
import glob
import csv
import multiprocessing as mp

def do_kernel_wrapper(args):
	return do_kernel(*args)

def do_kernel(k, tree1, tree2):
	return k.kernel(tree1, tree2)

if __name__ == '__main__':
	print("Initializing...")

	tree_dir = sys.argv[1]
	paths = glob.glob(tree_dir + "*.tre")
	tree_names = [s[len(tree_dir):-4] for s in paths]

	trees = [Phylo.read(x, 'newick') for x in paths]

	k2 = K2.PhyloKernel(scaleFactor=0.3, gaussFactor=2)

	for tree in trees:
		tree.ladderize()
		
	n = len(trees)
	
	p = mp.Pool()
	
	if sys.argv[2] != "":
		k = K2.PhyloKernel(scaleFactor=0.3, gaussFactor=2, labelFactor=1, labelFilter=sys.argv[2])
		print("Computing kernels (labelled)...")

		with open(tree_dir + "kernels.csv", 'w') as f:		
			w = csv.writer(f)
	
			w.writerow(tree_names)
		
			for i in xrange(0, n):
				print i + 1,
				sys.stdout.flush()
			
				kernels = p.map(do_kernel_wrapper, [(k, trees[i], trees[j]) for j in xrange(i, n)])
				
				row = [0] * i
				row.extend([str(x) for x in kernels])
				
				w.writerow(row)
				f.flush()
			
	print("\nComputing kernels (unlabelled)...")

	with open(tree_dir + "kernels.u.csv", 'w') as f:
		w = csv.writer(f)

		w.writerow(tree_names)
		
		for i in xrange(0, n):
			print i + 1,
			sys.stdout.flush()
			
			kernels = p.map(do_kernel_wrapper, [(k2, trees[i], trees[j]) for j in xrange(i, n)])
				
			row = [0] * i
			row.extend([str(x) for x in kernels])
	
			w.writerow(row)
			f.flush()