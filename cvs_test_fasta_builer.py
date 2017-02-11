#!/usr/bin/env python

import argparse
import sys
import numpy as np

possible = "ACDEFGHIKLMNPQRSTVWY"


def random_seq_from_dirichlet(bool_pow2, bool_even, bool_descend):
    if bool_pow2:
        alpha = []
        for a in range(0, 20):
            alpha.append(pow(2, 20) / pow(2.0, a))
    elif bool_even:
        alpha = np.full(20, 1)
        alpha[0] = 10

    for i in range(0, 1000):
        test_fasta = []
        for n in range(0, 20):

            if bool_descend:
                alpha = np.full(20, n + 1)
                alpha[0] = 100
                pmf = list(np.random.dirichlet(alpha, size=1)[0])
                test_fasta.append(np.random.choice(list(possible), p=pmf))
            else:
                pmf = list(np.random.dirichlet(alpha[:n + 1], size=1)[0])
                test_fasta.append(np.random.choice(list(possible[:n + 1]), p=pmf))
        fasta = ''.join(test_fasta)
        sys.stdout.write('{0}\n'.format(fasta))


def random_seq_from_uniform():
    for i in range(0, 1000):
        test_fasta = []
        for n in range(0, 20):
            test_fasta.append(np.random.choice(list(possible[:n + 1])))
        fasta = ''.join(test_fasta)
        sys.stdout.write('{0}\n'.format(fasta))


def build_test_fasta(bool_dirichlet_pow2, bool_dirichlet_even, bool_dirichlet_descend, seed):
    np.random.seed(seed)
    if bool_dirichlet_pow2 or bool_dirichlet_even or bool_dirichlet_descend:
        random_seq_from_dirichlet(bool_dirichlet_pow2, bool_dirichlet_even, bool_dirichlet_descend)
    else:
        random_seq_from_uniform()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""quick script to generate test cases for cvs analysis
        script generates 100, 20 amino acid long, fastas. first position is fully conserved. last is completely random.
        intermediary positions increase in randomness and amino acid identiy is picked uniformly by default"""
    )
    dirichlet_group = parser.add_mutually_exclusive_group()
    dirichlet_group.add_argument('-d_pow2', '--dirichlet_power2', action='store_true', default=False,
                                 help='pick amino acid for each position according to '
                                      'dirichlet distribution in stead of uniform. each amino acid considered has 1/2 '
                                      'less chance of being picked as the one before it')
    dirichlet_group.add_argument('-d_even', '--dirichlet_even', action='store_true', default=False,
                                 help='dirichlet distribution. first aa in possible heavily weighted. all others '
                                      'given approximately even weights')
    dirichlet_group.add_argument('-d_d', '--dirichlet_descend', action='store_true', default=False,
                                 help='dirichlet distribution. all 20 aas considered at every poisition. '
                                      'as postition number increase, the probability of picking the first aa in '
                                      'possible decreases, resulting in increasing noise')
    parser.add_argument('-s', '--seed', type=int, default=0,
                        help='set seed value used to generate random selections. same seed will give same selections')
    args = parser.parse_args()
    build_test_fasta(args.dirichlet_power2, args.dirichlet_even, args.dirichlet_descend, args.seed)
