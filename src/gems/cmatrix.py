# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

import collections
import itertools

import numpy as np

import gems.branched_matrix

def _block(n, m):
    a = np.eye(n, dtype=int)
    b = np.empty((n,m), dtype=int)
    q = m // n
    for i in range(q):
        b[:,i::q] = a
    return(b)


def _linear(n):
    return np.eye(n, k=1, dtype=int) + np.eye(n, k=-1, dtype=int)


def _pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def _group(iter1, iter2):
    if not len(iter2) % len(iter1):
        raise ValueError
    q = len(iter2) // len(iter1)
    it = iter(iter2)
    for i in iter1:
        yield (i, [next(it) for _ in range(q)])


def connectivity_matrix(expanded_symmetry):
    c = collections.Counter(expanded_symmetry)
    v = tuple(sorted(c.values(), reverse=True))

    # case I. macrocycle
    if len(c) == 1:
        letter, n = c.popitem()
        retm = _linear(n)
        retm[0,-1] = 1
        retm[-1,0] = 1

    # case II. linear molecule
    elif all(_ == 1 for _ in v):
        n = len(c)
        retm = _linear(n)

    # case III. branched
    elif all(i%j==0 for i,j in zip(v[:-1], v[1:])):
        retm = gems.branched_matrix.branched_molecule(expanded_symmetry)
        # v = tuple(sorted(c.values()))
        # a = (sum(v[:i]) for i in range(1+len(v)))
        # slices = {l: slice(i,j) for l, (i, j) in zip(slabels, _pairwise(a))}
        # if v[0] == 2:
        #     m[:2,:2] = np.eye(2, dtype=int)[:,::-1]

        # for a, b in _pairwise(slabels):
        #     sa = slices[a]
        #     sb = slices[b]
        #     block = _block(c[a],c[b])
        #     m[sa, sb] = block
        #     m[sb, sa] = block.T
        # retm = m
    else:
        print("molecule symmetry is not procesable")
        print(f"{count[b]} not divisible by {count[a]}")
        print("If this is a real molecule, you have to input the",
              "connectivity manually.")
        raise ValueError
    return retm
