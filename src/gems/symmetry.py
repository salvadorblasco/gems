# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2021 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""
Routines for isomorphism calculations.
"""

import collections
import itertools
import re

import numpy as np

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Molecule name identifier
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

def expand_name(input_name):
    """Expand the name of the molecule.

    Given a molecule symmetry in the short form, expand it
    to replace the numbers expressing multiplicity with that
    multiplicity.

    >>> expand_name('A2B4')
    'AABBB'
    """
    # import collections
    tcount = ''
    prev = None
    out = []
    for l in input_name:
        if l.isdigit():
            tcount += l
        else:
            if tcount:
                count = int(tcount)
                out.extend(int(tcount) * prev)
                tcount = ''
            else:
                if prev:
                    out.append(prev)
            prev = l

    if tcount:
        count = int(tcount)
        out.extend(count * prev)
    else:
        out.append(prev)
    out.sort()
    return "".join(out)


def expand_sort_name(input_name):
    expr = re.compile("(\w)(\d+)?") 
    matches = expr.finditer(input_name)
    mdict = {x[1].upper(): (1 if x[2] is None else int(x[2])) for x in matches}
    return "".join(mdict[a]*a for a in sorted(mdict, key=lambda x:mdict[x]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#    Molecular symmmetry
# ~~~~~~~~~~~~~~~~~~~~~~~~~~


def clasify_microstates(microstates, isomorphisms):
    """
    >>> isomorphisms = [array([[1, 0, 0],
    ...                        [0, 0, 1],
    ...                        [0, 1, 0]])]
    >>> microstates = {(1, 0, 0), (0, 0, 1), (0, 1, 0)}
    >>> clasify_microstates(microstates, isomorphisms)
    {((0, 1, 1),), ((1, 0, 1), (1, 1, 0))}
    """
    clasify = {ms: [ms] for ms in microstates}
    
    for ms, isomorphism in itertools.product(microstates, isomorphisms):
        current_list = clasify[ms]
        trr = tuple(np.dot(isomorphism, np.array(ms, dtype=np.int)))
        if trr not in current_list:
            current_list.append(trr)
    
    return {tuple(sorted(_)) for _ in clasify.values()}


def generate_microstates(molecule, level):
    """
    >>> generate_microstates("AAB", 1)
    {(1, 0, 0), (0, 0, 1), (0, 1, 0)}
    >>> generate_microstates("AAB", 2)
    {(1, 0, 1), (1, 1, 0), (0, 1, 1)}
    """
    size = len(molecule)
    return set(itertools.permutations(level*(1,)+(size-level)*(0,), size)) 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Isomorphisms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_isomorphisms(adjacency):
    trivial = _trivial_isomorphisms(adjacency)
    if trivial is None:
        size = len(adjacency)
        iterpermut = itertools.permutations(np.hsplit(np.eye(size, dtype=np.int), size), size)
        permutation_matrices = (np.hstack(c) for c in iterpermut)
        isomorphisms = []
        for pmatrix in permutation_matrices:
            trmatrix = np.dot(np.dot(pmatrix, adjacency), pmatrix.T)
            if np.all(trmatrix == adjacency) and not np.all(pmatrix == np.eye(size)):
                isomorphisms.append(pmatrix)
    else:
        isomorphisms = trivial
    return isomorphisms


def connectivity_matrix(molecule):
    size = len(molecule)

    if size < 3:
        m = np.ones((size, size), dtype=np.int)
    else:
        count = collections.Counter(molecule)
        slabels = list(sorted(count, key=lambda x: count[x]))
        m = np.zeros((size, size), dtype=np.int)

        start = 0
        slices = dict()
        for label in slabels:
            end = start + count[label]
            slices[label] = slice(start, end)
            start = end

        s = slices[slabels[0]]
        m[s, s] = 1
        for a, b in _pairwise(slabels):
            sa = slices[a]
            sb = slices[b]
            if count[b] % count[a] != 0:
                raise ValueError("molecule symmetry is not procesable"
                                 f"{count[b]} not divisible by {count[a]}"
                                 "If this is a real molecule, you have to input the"
                                 "connectivity manually.")
            mult = count[b]//count[a]
            sub = _subeye(count[a], mult)
            m[sa, sb] = sub
            m[sb, sa] = sub.T
    return m


def _pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _subeye(size, mult):
    assert isinstance (mult, int), mult
    if mult == 1:
        retv = np.eye(size, dtype=np.int)
    else:
        aux1 = np.hsplit(np.eye(size, dtype=np.int), size)
        aux2 = list(itertools.chain(tuple(itertools.repeat(b,mult)) for b in aux1)) 
        retv = np.hstack(list(itertools.chain.from_iterable(x for x in aux2)))
    return retv


def _trivial_isomorphisms(adjacency):
    if np.all(adjacency == np.ones_like(adjacency, dtype=np.int)):
        retv = len(adjacency) * (adjacency,)
    elif np.all(adjacency == np.eye(len(adjacency), dtype=np.int)):
        retv = (adjacency, )
    else:
        retv = None
    return retv
