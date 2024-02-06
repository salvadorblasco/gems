# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

import itertools
import re
import sys

import numpy as np

import gems.cmatrix

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _trivial_isomorphisms(adjacency):
    if np.all(adjacency == np.ones_like(adjacency, dtype=np.int)):
        # retv = len(adjacency) * (adjacency,)
        retv = len(adjacency) * (np.eye(len(adjacency), dtype=np.int),)
    elif np.all(adjacency == np.eye(len(adjacency), dtype=np.int)):
        retv = (adjacency, )
    else:
        retv = None
    return retv


def _permute_labels(labels, pmatrix):
    return "".join(labels[j] for j in pmatrix.nonzero()[-1])


def find_isomorphisms(adjacency, labels):
    """Given an adjacency matrix, find all isomorphisms by brute force.
    """
    def aux(pm):
        trlabels = _permute_labels(labels, pm)
        if trlabels != labels:
            return False, pm
        trmatrix = np.dot(np.dot(pm, adjacency), pm.T)
        return np.all(trmatrix == adjacency), pm

    size = len(adjacency)
    iterpermut = itertools.permutations(np.hsplit(np.eye(size, dtype=np.int), size), size)
    permutation_matrices = (np.hstack(c) for c in iterpermut)
    mymap = map(aux, permutation_matrices)
    isomorphisms = [m for c, m in mymap if c]
    return isomorphisms


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


# deprecated. use libuk.generate_microstates
def generate_microstates(m, level):
    """
    >>> generate_microstates(3, 1)
    {(1, 0, 0), (0, 0, 1), (0, 1, 0)}
    >>> generate_microstates(3, 2)
    {(1, 0, 1), (1, 1, 0), (0, 1, 1)}
    """
    return set(itertools.permutations(level*(1,)+(m-level)*(0,), m)) 


def subeye(size, mult):
    assert isinstance(mult, int), mult
    if mult == 1:
        retv = np.eye(size, dtype=np.int)
    else:
        aux1 = np.hsplit(np.eye(size, dtype=np.int), size)
        aux2 = list(itertools.chain(tuple(itertools.repeat(b, mult)) for b in aux1))
        retv = np.hstack(list(itertools.chain.from_iterable(x for x in aux2)))
    return retv


def expand_name(input_name):
    return re.sub(r'(\w)(\d+)', lambda x: x.group(1)*int(x.group(2)), input_name)


# # DELETE
# def repr_connectivity(connectivity, molecule):
#     print("  ", "".join(f"{m:2}" for m in molecule))
#     for rlet, row in zip(molecule, connectivity):
#         print(f"{rlet}", "".join(f"{x:2}" for x in row))
# 
# def repr_clasified_microstates(my_set):
#     print('{', end='')
#     for t in my_set:
#         print('(', end='')
#         print(",".join(('(' + "".join(str(_) for _ in tt) + ')' for tt in t)), end='')
#         print(')', end=',')
#     print('}')
# 
# 
# def ms_to_text(molecule, ms):
#     return "".join(molecule[n].upper() if b else molecule[n].lower() for n, b in enumerate(ms))
# 
# 
# def print_isomorphisms(isomorphisms, molecule):
#     for n, isomorphism in enumerate(isomorphisms):
#         print(f'isomorphism #{n}')
#         repr_connectivity(isomorphism, molecule)
# 
# 
# # DELETE
# def print_molecule_report(input_molecule):
#     print(f"Input molecule: {input_molecule}")
#     print()
#     molecule = expand_name(input_molecule)
#     print("Expanded molecule: ", molecule)
#     matrix1 = cmatrix.connectivity_matrix(molecule)
#     print("Tentative connectivity matrix")
#     repr_connectivity(matrix1, molecule)
#     print()
#     isomorphisms = find_isomorphisms(matrix1, molecule)
#     print(f"found {len(isomorphisms)} isomorphisms")
#     print()
#     # print_isomorphisms(isomorphisms, molecule)
#     # print()
# 
#     total_microstates = 0
#     parameters = []
#     for level in range(1+len(molecule)):
#         print(5*'-', level, 5*'-')
#         microstates = generate_microstates(len(molecule), level)
#         # microstates = generate_microstates(molecule, level)
#         my_set = clasify_microstates(microstates, isomorphisms)
#         repr_clasified_microstates(my_set)
#         total_microstates += len(my_set)
#         if 0 < level < 4:
#             parameters.append(len(my_set))
#         print(f"{len(my_set)} microstates")
#         print("m: equivalent microstates")
#         aux = [(len(group),
#                 (ms_to_text(molecule, ms) for ms in group))
#                for group in my_set]
#         for a, b in sorted(aux, key=lambda x: x[0]):
#             print(f"{a}:", ", ".join(b))
#         print()
# 
#     print(f"Total of {total_microstates} independent microstates")
#     for tx, nn in zip(('first', 'second', 'third'), parameters):
#         print(f"Total {tx}-level parameters: {nn}")
