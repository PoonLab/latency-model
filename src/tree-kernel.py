import argparse
from Bio import Phylo
import sys
import math
import re

def preprocess(tree, normalize='mean', ladderize=True, labelRegex=None):
    if ladderize:
        tree.ladderize()
    normalize_tree(tree, normalize)
    annotate_tree(tree, labelRegex)

def normalize_tree(t, mode='mean'):
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
        branch_lengths = [branch.branch_length for branch in
                          branches[int(not t.rooted):]]
        branch_lengths.sort()

        if nbranches % 2 == 0:
            median_branch_length = (branch_lengths[(nbranches / 2) - 1] +
                                    branch_lengths[nbranches / 2]) / 2.
        else:
            median_branch_length = branch_lengths[nbranches / 2]

        for branch in branches[int(not t.rooted):]:
            branch.branch_length /= median_branch_length
    elif mode == 'none':
        pass
    else:
        print 'ERROR: Unrecognized mode in normalize_tree():', mode
        sys.exit()


0
def annotate_tree(t, labelRegex=None):
    """ Update Clade objects with productions """
    for tip in t.get_terminals():
        tip.production = 0
        if labelRegex:
            matches = labelRegex.findall(tip.name)
            tip.label = matches[0] if matches else ''

    for i, node in enumerate(t.get_nonterminals(order='postorder')):
        children = node.clades
        nterms = sum([c.production == 0 for c in children])
        node.production = nterms + 1
        node.index = i
        branch_lengths = [c.branch_length for c in node.clades]
        node.bl = branch_lengths
        node.sqbl = sum([bl ** 2 for bl in branch_lengths])


def kernel(t1, t2, decayFactor, gaussFactor, sigma, labelFilter, labelFactor):
    settings = {
        'sigma': sigma,
        'decayFactor': decayFactor,
        'gaussFactor': gaussFactor,
        'labelFactor': labelFactor,
        'labelFilter': labelFilter,
    }
    dp_matrix = {}  # dynamic programming
    k = 0.  # result

    nodes1 = t1.get_nonterminals(order='postorder')
    if not hasattr(nodes1[0], 'production'):
        annotate_tree(t1)
    nodes2 = t2.get_nonterminals(order='postorder')
    if not hasattr(nodes2[0], 'production'):
        annotate_tree(t2)

    for ni, n1 in enumerate(nodes1):
        for n2 in nodes2:
            if n1.production == n2.production:
                if settings['labelFilter']:
                    result = labeled_kernel(n1, n2, settings, dp_matrix)
                else:
                    result = unlabeled_kernel(n1, n2, settings, dp_matrix)
                dp_matrix[(n1.index, n2.index)] = result
                k += result
    return k


def labeled_kernel(n1, n2, settings, dp_matrix):
    res1 = math.exp(-1. / settings['gaussFactor'] * (n1.sqbl + n2.sqbl - 2 * sum(
        [(n1.bl[i] * n2.bl[i]) for i in range(len(n1.bl))])))
    res2 = math.exp(-1. / settings['gaussFactor'] * (n1.sqbl + n2.sqbl - 2 * sum(
        [(n1.bl[i] * n2.bl[(i + 1) % 2]) for i in range(len(n1.bl))])))

    for cn1 in range(2):
        c1 = n1.clades[cn1]
        c21 = n2.clades[cn1]
        c22 = n2.clades[(cn1 + 1) % 2]

        if c1.production != c21.production:
            res1 *= settings['sigma']
        elif c1.production == 0:
            # branches are terminal
            labelWeight = 1 if c1.label == c21.label else settings['labelFactor']

            res1 *= settings['sigma'] + settings['decayFactor'] * labelWeight
        else:
            try:
                res1 *= settings['sigma'] + dp_matrix[(c1.index, c21.index)]
            except KeyError:
                res1 *= settings['sigma']

        if c1.production != c22.production:
            res2 *= settings['sigma']
        elif c1.production == 0:
            # branches are terminal
            labelWeight = 1 if c1.label == c22.label else settings['labelFactor']

            res2 *= settings['sigma'] + settings['decayFactor'] * labelWeight
        else:
            try:
                res2 *= settings['sigma'] + dp_matrix[(c1.index, c22.index)]
            except KeyError:
                res2 *= settings['sigma']

    return settings['decayFactor'] * max(res1, res2)


def unlabeled_kernel(n1, n2, settings, dp_matrix):
    res = settings['decayFactor'] * math.exp( -1./settings['gaussFactor'] * (n1.sqbl + n2.sqbl - 2*sum(
        [(n1.bl[i] * n2.bl[i]) for i in range(len(n1.bl))]
    )))
    for cn1 in range(2):
        c1 = n1.clades[cn1]
        c2 = n2.clades[cn1]
        if c1.production != c2.production:
            res *= settings['sigma']
        elif c1.production == 0:
            # branches are terminal (c2 == c1)
            res *= settings['sigma'] + settings['decayFactor']
        else:
            try:
                res *= settings['sigma'] + dp_matrix[(c1.index, c2.index)]
            except KeyError:
                res *= settings['sigma']
    return res


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculates subset tree kernel on phylogenies using radial '
                    'basis function (RBF) kernel to penalize differences in '
                    'branch lengths.'
    )
    parser.add_argument('t1', type=argparse.FileType('rU'),
                        help='<input> File containing Newick tree string(s)')
    parser.add_argument('t2', nargs='?', type=argparse.FileType('rU'), default=None,
                        help='<optional, input> Second file containing Newick tree '
                             'string for pairwise comparison.')
    parser.add_argument('out', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='<optional, output> File to write results to (default: stdout)')

    parser.add_argument('--decay', type=float, default=0.2,
                        help='Decay factor penalizing large subset trees.')
    parser.add_argument('--rbfSigma', type=float, default=1.0,
                        help='Variance parameter of RBF kernel.')
    parser.add_argument('--sigma', type=float, default=1.0,
                        help='sigma=1.0 -> subset tree kernel; sigma=0.0 -> subtree kernel')
    parser.add_argument('--labelFactor', type=float, default=1.0,
                        help='Factor for penalizing mismatched tip labels.')
    parser.add_argument('--labelFilter', type=str, default=None,
                        help='')

    parser.add_argument('--noLadderize', action='store_true',
                        help='If set, will not rotate branches so that the '
                             'branches with the most tips are always on the '
                             'same side.')
    parser.add_argument('--normalize', type=str, default='mean',
                        help='Branch length normalization mode (mean, median, none).')

    return parser.parse_args()



def main():
    args = parse_args()

    labelRegex = re.compile(args.labelFilter) if args.labelFilter else None
    ladderize = False if args.noLadderize else True

    if args.t2:
        # calculate kernel score for two trees only
        t1 = preprocess(Phylo.read(args.t1, 'newick'), args.normalize, ladderize, labelRegex)
        t2 = preprocess(Phylo.read(args.t2, 'newick'), args.normalize, ladderize, labelRegex)
        k = kernel(t1, t2, args.decay, args.rbfSigma, args.sigma, args.labelFilter, args.labelFactor)
        args.out.write(k)
    else:
        # calculate kernel matrix for all trees in first file
        trees = []
        for tree in Phylo.parse(args.t1, 'newick'):
            trees.append(preprocess(tree))


if __name__ == '__main__':
    main()
