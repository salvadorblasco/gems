#!/usr/bin/python3

import collections
import itertools

import numpy as np


def _pairwise2(iterable2):
    def f():
        yield None
        yield from iterable2

    def g():
        yield from iterable2
        yield None

    return zip(f(), g())


def _pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _group(iter1, iter2):
    if len(iter2) % len(iter1) > 0:
        raise ValueError
    q = len(iter2) // len(iter1)
    it = iter(iter2)
    for i in iter1:
        yield (i, [next(it) for _ in range(q)])


def branched_molecule(molecule):
    c = collections.Counter(molecule)
    kk = {ll:tuple(n for n,l in enumerate(molecule) if l==ll) for ll in set(molecule)}
    size = len(molecule)
    m = np.zeros((size, size), dtype=np.int)

    for a, b in _pairwise2(sorted(c, key=c.get)):
        if a is None:
            ib = kk[b]
            if len(ib) == 2:
                m[ib[0], ib[1]] = 1
                m[ib[1], ib[0]] = 1
        elif b is None:
            break
        else:
            q = len(kk[b]) // len(kk[a])
            for aa, bb in _group(kk[a], kk[b]):
                m[aa, bb] = 1
                m[bb, aa] = 1
    return m
