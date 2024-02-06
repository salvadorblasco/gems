# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

""" Variable Mapping Cascade module """

import itertools
import sys

import gems.libuk


class FittingParams():
    def __init__(self, mapping):
        self.constraints = {}
        self.restraints = []
        self.restrained_ids = []
        self.restraint_values = []
        self.levels = {}
        self.size = len(tuple(mapping.keys())[0])
        self.mapping = mapping

        parm_keys = itertools.chain.from_iterable([tuple(gems.libuk.filter_by_macro(mapping, n)) for n in (1,2,3)])
        self.parameters = dict.fromkeys(set(parm_keys), 0.0)
        self.parameter_errors = dict.fromkeys(self.parameters.keys(), 0.0)

        self.__prepareids()

    def __prepareids(self):
        aux = {(sum(m[0]), min(sum(a*2**n for n, a in enumerate(b)) for b in m)): m
               for m in set(self.mapping.values())}
        self.ids = {n:aux[k] for n, k in enumerate(sorted(aux))}
        self.revid = {aux[k]:n for n, k in enumerate(sorted(aux))}

    def create_constraint(self, msid):
        microstate = self.ids[msid]
        self.constraints[microstate] = self.parameters.pop(microstate)

    def create_restraint(self, ids):
        aux = [self.parameters[self.ids[i]] for i in ids]
        values = [v/aux[0] for v in aux]
        self.restraints.append(ids)
        self.restrained_ids.extend(ids)
        self.restraint_values.append(values)
        for i in ids[1:]:
            del self.parameters[self.ids[i]]

    def error_free_energy(self, microstate, max_level=3):
        mult = (1, 1/2, 1/6)
        return self.__common_free_energy(self.evmc, mult, microstate, max_level) 

    def free_energy(self, microstate, max_level=3):
        mult = (-1, 1/2, 1/6)
        return self.__common_free_energy(self.vmc, mult, microstate, max_level) 

    def __common_free_energy(self, fvmc: callable, mult: tuple, microstate, max_level=3):
        energy = 0.0
        LN10 = 2.3025851
        for level in range(max_level):
            for ms in gems.isomorph.generate_microstates(self.size, 1+level):   # TODO change this function
                if sum(i*j for i,j in zip(ms, microstate)) == 1+level:
                    p = fvmc(ms)
                    energy += p * mult[level]
        return energy*LN10

    def name_microstate(self, microstate):
        return "".join(sorted(l for l, i in zip(self.molecule, microstate) if i))

    def microstate_id(self, microstate):
        c = self.mapping[microstate]
        return self.revid[c]

    def microstate_status(self, microstate):
        c = self.mapping[microstate]
        i = self.revid[c]
        if c in self.constraints:
            return "fixed"
        elif i in self.restrained_ids:
            # for n, r in enumerate(self.restraints):
            #     if i in r:
            #         break
            n, _, _ = self.__locate_restraint(c)
            return f"restrained[{n}]"
        else:
            return "refine"

    def evmc(self, microstate):
        return self.__common_vmc(microstate, True)

    def vmc(self, microstate):
        return self.__common_vmc(microstate, False)

    def __common_vmc(self, microstate, errors=False):
        if errors:
            parameters = self.parameter_errors
        else:
            parameters = self.parameters

        ms = self.mapping[microstate]
        if ms in parameters:
            return parameters[ms]
        elif ms in self.constraints:
            return 0.0 if errors else self.constraints[ms]
        elif (i := self.revid[ms]) in self.restrained_ids:
            # for n, r in enumerate(self.restraints):
            #     if i in r:
            #         j = r.index(i)
            #         break
            # ref = self.restraints[n][0]
            # value = parameters[self.ids[ref]] * self.restraint_values[n][j]
            _, _, value = self.__locate_restraint(microstate)
            return value
        else:
            raise ValueError(f'microstate {microstate} non-existent')

    def __locate_restraint(self, microstate):
        i = self.revid[microstate]
        for n, r in enumerate(self.restraints):
            if i in r:
                j = r.index(i)
                restr_n = n
                break
        ref = self.restraints[restr_n][0]
        value = self.parameters[self.ids[ref]] * self.restraint_values[restr_n][j]
        return restr_n, ref, value

    def set_values(self, values):
        self.__common_set(values, self.parameters)

    def set_errors(self, values):
        self.__common_set(values, self.parameter_errors)

    def __common_set(self, values, store):
        ival = itertools.chain(values, itertools.repeat(0.0))
        for _id in sorted(self.ids):
            microstate = self.ids[_id]
            if microstate in store:
                store[microstate] = next(ival)
