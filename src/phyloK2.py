# by Art FY Poon, 2012
# modified by Rosemary McCloskey, 2014
# modified by Bradley R Jones, 2016

from Bio import Phylo
from numpy import zeros
import math
import multiprocessing as mp
import re
import time

class PhyloKernel:
    def __init__(self,
                 kmat=None,
                 rotate='ladder',
                 rotate2='none',
                 subtree=False,
                 normalize='mean',
                 sigma=1,
                 gaussFactor=1,
                 withLengths=True,
                 decayFactor=0.1,
                 verbose=False,
                 resolve_poly=False,
                 labelFactor=1,
                 labelFilter=None,
                 **kwargs
        ):
        """
        requires a list of Phylo.Tree objects
        can cast iterator returned by Phylo.parse() as list
        """
        self.resolve_poly = resolve_poly
        self.normalize = normalize
        
        self.trees = []
        
        self.kmat = []
        self.is_kmat_computed = False
        
        # using **kwargs would probably make this cleaner
        self.rotate = rotate
        self.rotate2 = rotate2
        self.subtree = subtree
        self.normalize = normalize
        self.sigma = sigma
        self.gaussFactor = gaussFactor
        self.withLengths = withLengths
        self.decayFactor = decayFactor
        self.verbose = verbose
        self.resolve_poly = resolve_poly
        
        # labelling
        self.labelFactor = labelFactor
        self.labelFilter = labelFilter
        if labelFilter is not None:
            self.labelRegex = re.compile(self.labelFilter)
        
        self.pcache = {}
        self.subtrees = {} # used for matching polytomies
        #self.cache_productions()
        
        self.sigma = sigma
        self.gaussFactor = gaussFactor
        self.decayFactor = decayFactor
        self.withLengths = withLengths
        
        self.verbose = verbose
        
        if self.verbose:
            print('creating PhyloKernel with settings')
            print('sigma = {}'.format(self.sigma))
            print('gaussFactor = {}'.format(self.gaussFactor))
            print('decayFactor = {}'.format(self.decayFactor))
            print('labelFactor = {}'.format(self.labelFactor))
        
    @property
    def ntrees (self):
        return len(self.trees)

    def reset_matrix(self):
        self.kmat = zeros((self.ntrees, self.ntrees))

    def load_trees_from_file (self, handle):
        """
        Parse a file containing Newick tree strings
        """
        self.trees = []
        
        tree_iter = Phylo.parse(handle, 'newick')
        
        for t in tree_iter:
            if self.rotate=='ladder':
                t.ladderize()

            if self.normalize != 'none': self.normalize_tree(t, mode=self.normalize)

            self.annotate_tree(t)
            self.trees.append(t)
            
        #self.kmat = [[0 for i in self.ntrees] for j in self.ntrees]
        self.is_kmat_computed = False
        self.delta_values = {}
    
    def normalize_tree (self, t, mode='median'):
        """
        Normalize branch lengths in tree by mean branch length.
        This helps us compare trees of different overall size.
        Ignore the root as its branch length is meaningless.
        """
        
        # compute number of branches in tree
        branches = t.get_nonterminals() + t.get_terminals()
        nbranches = len(branches) - 1

        if mode == 'mean':
            if t.root.branch_length is None:
                t.root.branch_length = 0.
            tree_length = t.total_branch_length() - t.root.branch_length
            mean_branch_length = tree_length / nbranches
            
            for branch in branches[int(not t.rooted):]:
                branch.branch_length /= mean_branch_length
            
        elif mode == 'median':
            branch_lengths = [branch.branch_length for branch in branches[int(not t.rooted):]]
            branch_lengths.sort()
            
            if nbranches%2 == 0:
                median_branch_length = (branch_lengths[(nbranches/2)-1] + 
                                        branch_lengths[nbranches/2]) / 2.
            else:
                median_branch_length = branch_lengths[nbranches/2]
            
            for branch in branches[int(not t.rooted):]:
                branch.branch_length /= median_branch_length

    def annotate_tree(self, t):
        """
        Add annotations to Clade objects in place
        """
        for tip in t.get_terminals():
            tip.production = 0
            
            if self.labelFilter != None:
                matches = self.labelRegex.findall(tip.name)
                tip.label = matches[0] if matches else ''

        for i, node in enumerate(t.get_nonterminals(order='postorder')):
            children = node.clades
            nterms = sum( [c.production == 0 for c in children] )
            node.production = nterms + 1
            node.index = i
            branch_lengths = [c.branch_length for c in node.clades]
            node.bl = branch_lengths
            node.sqbl = sum([bl**2 for bl in branch_lengths])
            	

    def compute_matrix(self):
        for i in range(self.ntrees):
            for j in range(i, self.ntrees):
                self.kmat[i,j] = self.kmat[j,i] = self.kernel(self.trees[i], self.trees[j])
                    
                if self.verbose:
                    print('%d\t%d\t%f' % (i, j, self.kmat[i,j]))
        
        self.is_kmat_computed = True

#	def parse_label(self, clade):
		

    def kernel(self, t1, t2, myrank=None, nprocs=None, dp_matrix={}):
        """
        Recursive function for computing tree convolution
        kernel.  Adapted from Moschitti (2006) Making tree kernels
        practical for natural language learning. Proceedings of the 
        11th Conference of the European Chapter of the Association 
        for Computational Linguistics.
        """
        nodes1 = t1.get_nonterminals(order='postorder')
        nodes2 = t2.get_nonterminals(order='postorder')
        k = 0
        
        if not hasattr(nodes1[0], 'production'):
            self.annotate_tree(t1)
        if not hasattr(nodes2[0], 'production'):
            self.annotate_tree(t2)

        # iterate over non-terminals, visiting children before parents
        for ni, n1 in enumerate(nodes1):
            if myrank is not None and nprocs is not None and ni % nprocs != myrank:
                continue
                        
            # if multithreading wait for children to be processed first
            if myrank is not None:
#                print("in: " + str(n1.index))
                
                for c in n1.clades:
                    if hasattr(c, 'index'):
                        while c.index not in dp_matrix:
#	    	        		print(str(n1.index) + " wait on: " + str(c.index))
    		        		time.sleep(0.0001)

            for n2 in nodes2:
                if n1.production == n2.production:
                    if self.labelFilter is not None:                    
                        res1 = math.exp( -1. / self.gaussFactor * (n1.sqbl + n2.sqbl - 2*sum([(n1.bl[i]*n2.bl[i]) for i in range(len(n1.bl))]))) 
                        res2 = math.exp( -1. / self.gaussFactor * (n1.sqbl + n2.sqbl - 2*sum([(n1.bl[i]*n2.bl[(i+1)%2]) for i in range(len(n1.bl))])))
                        
                        for cn1 in range(2):
                            c1 = n1.clades[cn1]
                            c21 = n2.clades[cn1]
                            c22 = n2.clades[(cn1 + 1) % 2]
                        
                            if c1.production != c21.production:
                                res1 *= self.sigma
                            elif c1.production == 0:
                                # branches are terminal
                                labelWeight = 1 if c1.label == c21.label else self.labelFactor
                            	
                                res1 *= self.sigma + self.decayFactor * labelWeight
                            else:
                                try:
                                    res1 *= self.sigma + dp_matrix[(c1.index, c21.index)]
                                except KeyError:
                                    res1 *= self.sigma
                                
                            if c1.production != c22.production:
                                res2 *= self.sigma
                            elif c1.production == 0:
                                # branches are terminal
                                labelWeight = 1 if c1.label == c22.label else self.labelFactor
                            	
                                res2 *= self.sigma + self.decayFactor * labelWeight
                            else:
                                try:
                                    res2 *= self.sigma + dp_matrix[(c1.index, c22.index)]
                                except KeyError:
                                    res2 *= self.sigma
                                                        
                        res = self.decayFactor * max(res1, res2)
                    else:
                        res = self.decayFactor * math.exp( -1. / self.gaussFactor * (n1.sqbl + n2.sqbl - 2*sum([(n1.bl[i]*n2.bl[i]) for i in range(len(n1.bl))]))) 
                    
                        for cn1 in range(2):
                            c1 = n1.clades[cn1]
                            c2 = n2.clades[cn1]
                        
                            if c1.production != c2.production:
                                res *= self.sigma
                            elif c1.production == 0:
                                # branches are terminal
                                res *= self.sigma + self.decayFactor
                            else:
                                try:
                                    res *= self.sigma + dp_matrix[(c1.index, c2.index)]
                                except KeyError:
                                    res *= self.sigma
            
                    dp_matrix[(n1.index, n2.index)] = res
                    k += res
                
            if myrank is not None:    
                dp_matrix[n1.index] = 0
            
#            if myrank is not None:
#                print("out: " + str(n1.index))

#        if output is None:
        return k
        
            
    def kernel_parallel(self, t1, t2, nthreads):
        """
        Wrapper around kernel().
        Attempt to use Python multiprocessing module to speed up computation.
        Borrowing code snippets from:
          http://sebastianraschka.com/Articles/2014_multiprocessing_intro.html

        :param t1: first Phylo.Tree to be compared
        :param t2: second Phylo.Tree to be compared
        :param nthreads: number of threads in pool
        :return: kernel score (double)
        """
                
        manager = mp.Manager()
        dp_matrix = manager.dict()
        processes = [mp.Process(target=self.kernel, args=(t1, t2, i, nthreads, dp_matrix))
                     for i in range(nthreads)]
        for p in processes:
            p.start()

        # exit completed processes
        for p in processes:
            p.join()

        # collect results and calculate sum
        k = sum(dp_matrix.values())
        return k
