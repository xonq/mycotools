#! /usr/bin/env python3

"""
This is a pairwise subsampling procedure adapted from the following
publication:

Emmanuel Sapin, Matthew C Keller, Novel approach for parallelizing pairwise
comparison problems as applied to detecting segments identical by decent in
whole-genome data, Bioinformatics, Volume 37, Issue 15, 1 August 2021, Pages
2121–2125, https://doi.org/10.1093/bioinformatics/btab084 

Script authored by Zachary Konkel
"""

import sys
from math import ceil
from itertools import combinations

def isprime(i):
    for i0 in range(2, ceil(i/2) + 1): # from 2 to ceiling of i / 2
        if i % i0 == 0: # if it is completely divisible 
            return False # it is not prime
    return True # it is prime

def lowest_square_of_prime(i):
    for i0 in range(2, ceil(i)):
        if isprime(i0):
            if i0**2 > i:
                return i0
    return False

def fill_groups(n, ntop, l_sq_p):
    group, comb_data, out = [], [], []
    for i0 in range(0, ntop, l_sq_p): #
        for i1 in range(i0, i0 + l_sq_p):
            if i1 <= n:
                group.append(i1)
        comb_data.extend([x for x in combinations(sorted(group), 2)])
        out.append(group)
        group = []

    for i0 in range(0, l_sq_p):
        for i1 in range(0, l_sq_p):
            check = i0 + i1 * l_sq_p
            if check < n:
                group.append(check)
        comb_data.extend([x for x in combinations(sorted(group), 2)])
        out.append(group)
        group = []

    moveup = [0 for i in range(l_sq_p)]
    for increment in range(1, l_sq_p):
        for i0 in range(1, l_sq_p):
            moveup[i0] = (i0 * increment) % l_sq_p
        for i0 in range(l_sq_p):
            for i1 in range(l_sq_p):
                check = (i0 + moveup[i1]) % l_sq_p + i1 * l_sq_p
                if check <= n:
                    group.append(check)
            comb_data.extend([x for x in combinations(sorted(group), 2)])
            out.append(group)
            group = []

    return comb_data, out

def final_group(comb_data, n, groups):
    shouldbe = set([x for x in combinations(range(n + 1), 2)])
    observed = set(comb_data)
    missing = list(shouldbe.difference(observed))

    missing_i = []
    for miss_comb in missing:
        [missing_i.append(x) for x in miss_comb]

    groups.append([x for x in list(set(missing_i))])
    return groups

def main(n):
    """n = int(number of samples)"""
    l_sq_p = lowest_square_of_prime(n)
    # number of persons per group
    ntop = l_sq_p * l_sq_p # lowest squared prime > samples
    # number of elements in the matrix
    nb_group_per_person = l_sq_p + 1 # number of group per person
    nb_group = ntop + l_sq_p # number of groups
    comb_data, out = fill_groups(n, ntop, l_sq_p)
    out = final_group(comb_data, n, out)
    return [sorted(x) for x in out if x]

def cli():
    usage = 's2subs.py <SAMPLE_NUMBER> [OUTPUT_FILE]'
    if not {'-h', '--help', '-help'}.isdisjoint(set(sys.argv)) \
        or len(sys.argv) < 2:
        print('\n' + usage + '\n')
        sys.exit(1)

    n = int(sys.argv[1])
    out = main(n)

    if len(sys.argv) == 2:
        print('\n'.join([
            ' '.join([str(y) for y in x])
            for x in out
            ]))
    else: # there is an output file
        with open(sys.argv[2], 'w') as handle:
            handle.write('\n'.join([
                ' '.join([str(y) for y in x])
                for x in out
                ]))
    sys.exit(0)


if __name__ == '__main__':
    cli()
