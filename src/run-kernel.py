"""
Use MPI to rapidly compute kernel matrix on simulated trees

test.kmat.csv : no labelFilter (all tips have same labels)
test.lab-kmat.csv : labelFilter='\.(\w+)$', tips are labelled by deme; labelFactor=1
                    I think this means that mismatches are not penalized...
test.lab-kmat0.csv : labelFilter='\.(\w+)$', tips are labelled by deme; labelFactor=0
"""

from phyloK2 import PhyloKernel as PK
from glob import glob
from Bio import Phylo
import math
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
nprocs = comm.Get_size()

k = PK(decayFactor=0.2, gaussFactor=2, verbose=True, labelFactor=0, labelFilter='\.(\w+)$')
#k = PK(decayFactor=0.2, gaussFactor=2, verbose=True)

files = glob('../data/test/*.nwk')


# gather trees from different files
print 'loading trees'
treatments = []
for file in files:
    aL = float(file.split('/')[-1].strip('aL-.nwk'))
    #counter = 0
    for tree in Phylo.parse(file, 'newick'):
        tree.root.branch_length = 0.
        tree.ladderize()
        k.normalize_tree(tree, 'mean')
        k.annotate_tree(tree)
        k.trees.append(tree)
        treatments.append(aL)
        #counter += 1
        #if counter == 10:
        #    break

print treatments
sys.exit()
# initialize kernel matrix

k.reset_matrix()
#k.compute_matrix()
kdict = {}
for i in range(k.ntrees):
    for j in range(i, k.ntrees):
        index = i*k.ntrees + j
        if index % nprocs != my_rank:
            continue

        score = k.kernel(k.trees[i], k.trees[j])
        kdict.update({(i,j): score})
        print '(', my_rank, 'of', nprocs, ')', i, j, score

comm.Barrier()
kdicts = comm.gather(kdict, root=0)

if my_rank == 0:
    final = {}
    for kdict in kdicts:
        final.update(kdict)

    # convert into square matrix
    for i in range(k.ntrees):
        for j in range(i, k.ntrees):
            k.kmat[i,j] = k.kmat[j,i] = final[(i,j)]

    # normalize kernel scores and output
    with open('test.kmat0.csv', 'w') as outfile:
        for i in range(k.ntrees):
            row = [k.kmat[i,j] / math.sqrt(k.kmat[i,i]*k.kmat[j,j]) for j in range(k.ntrees)]
            outfile.write(','.join(map(str, row)))
            outfile.write('\n')
