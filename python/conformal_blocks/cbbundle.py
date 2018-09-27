from __future__ import division
import math, fractions, itertools, functools
from conformal_blocks.lie import SimpleLieAlgebra, TypeALieAlgebra, TypeBLieAlgebra, TypeCLieAlgebra, TypeDLieAlgebra, _Root
try:
    import sage.all as sage
    def Fraction(x,y):
        try:
            return sage.Rational((x, y))
        except TypeError:
            return x/y
except ImportError:
    from fractions import Fraction
'''
Created on Nov 10, 2016

@author: mjschust
'''

class ConformalBlocksBundle(object):
    """
    A class representing a conformal blocks vector bundle.
    """

    def __init__(self, liealg, weights, level):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param weights: A list of tuples of integers: the weights of the conformal blocks bundle.
        :param level: A positive integer: the level of the conformal blocks bundle.
        """
        self.liealg = liealg
        new_weights = []
        for wt in weights:
            new_weights.append(tuple(wt))
        self.weights = new_weights
        self.level = level
        self._rank = -1

    def get_rank(self):
        """
        Computes the rank of the conformal blocks bundle.  The algorithm uses factorization, then
        the fusion product to compute the 3-point ranks.

        :return: An integer: the rank of the bundle.
        """
        if self._rank < 0:
            self._rank = self._get_rank(self.weights, self.level)

        return self._rank

    def _get_rank(self, weights, level):
        """
        Computes the rank of the conformal blocks bundle with given weights and level.
        The algorithm uses the fusion product and factorization.

        :param weights: A list of tuples of integers: the list of weights.
        :param level: A positive integer: corresponds to the level of the fusion product.
        :return: An integer: the rank of the bundle.
        """
        # Find weights with largest and smallest corresponding rep's
        liealg = self.liealg
        min_dim = max_dim = liealg.get_rep_dim(weights[0])
        min_index = max_index = 0
        for i in range(len(weights)):
            dim = liealg.get_rep_dim(weights[i])
            if dim < min_dim:
                min_dim = dim
                min_index = i
            if dim > max_dim:
                max_dim = dim
                max_index = i
        # Covers the case when all dimensions are the same
        if min_index == max_index:
            max_index = min_index + 1

        fus_prod = liealg.fusion(weights[min_index], weights[max_index], level)
        # indices = min_index, max_index
        # factor_list = [wt for (i, wt) in enumerate(weights) if i not in indices]
        factor_list = []
        for i in range(len(weights)):
            if i != min_index and i != max_index:
                factor_list.append(weights[i])
        multi_fus_prod = liealg.multi_fusion(factor_list, level)

        ret_val = 0
        for mu_star in fus_prod:
            mult = fus_prod[mu_star]
            mu = liealg.get_dual_weight(mu_star)
            if mu in multi_fus_prod:
                ret_val += mult * multi_fus_prod[mu]

        return ret_val

    #Original version of the above method.  Uses less memory but runs an order of magnitude slower.
    def _alt_compute_rank(self, weights, level):
        # Find weights with largest and smallest corresponding rep's
        liealg = self.liealg
        min_dim = max_dim = liealg.get_rep_dim(weights[0])
        min_index = max_index = 0
        for i in range(len(weights)):
            dim = liealg.get_rep_dim(weights[i])
            if dim < min_dim:
                min_dim = dim
                min_index = i
            if dim > max_dim:
                max_dim = dim
                max_index = i
        # Covers the case when all dimensions are the same
        if min_index == max_index:
            max_index = min_index + 1

        fus_prod = liealg.fusion(weights[min_index], weights[max_index], level)
        # indices = min_index, max_index
        # factor_list = [wt for (i, wt) in enumerate(weights) if i not in indices]
        factor_list = []
        for i in range(len(weights)):
            if i != min_index and i != max_index:
                factor_list.append(weights[i])

        # Three point case is given by the fusion product
        if len(factor_list) == 1:
            dual_wt3 = liealg.get_dual_weight(factor_list[0])
            if dual_wt3 in fus_prod:
                return fus_prod[dual_wt3]
            else:
                return 0

        # If more than three points, factor
        ret_val = 0
        for wt in fus_prod:
            mult = fus_prod[wt]
            if mult > 0:
                ret_val = ret_val + mult * self._alt_compute_rank(factor_list + [wt], level)

        return ret_val

    def get_symmetrized_divisor(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle.

        :return: A list of numbers: the divisor given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """

        ret_val = []
        n = len(self.weights)
        weighted_rank = 0
        for wt in self.weights:
            weighted_rank += self.liealg.casimirScalar(wt)
        if self.liealg.exact:
            weighted_rank = Fraction(self.get_rank() * weighted_rank, n * (n - 1))
        else:
            weighted_rank = self.get_rank() * weighted_rank / (n * (n - 1))

        point_indices = [i for i in range(0, n)]
        for i in range(2, n // 2 + 1):
            coord = i * (n - i) * weighted_rank
            sum = 0
            for subset in itertools.combinations(point_indices, i):
                #Could be more efficient here
                wt_list1 = []
                wt_list2 = []
                for j in range(0,n):
                    if j in subset:
                        wt_list1.append(self.weights[j])
                    else:
                        wt_list2.append(self.weights[j])

                prod = self.liealg.multi_fusion(wt_list1, self.level)
                for mu_star in prod.keys():
                    mu = self.liealg.get_dual_weight(mu_star)
                    sum += self.liealg.casimirScalar(mu) * self._get_rank(wt_list1 + [mu], self.level) * self._get_rank(wt_list2 + [mu_star], self.level)

            if self.liealg.exact:
                sum = Fraction(sum * math.factorial(i) * math.factorial(n - i), math.factorial(n))
                coord = Fraction(coord - sum, 2 * (self.level + self.liealg.dual_coxeter()))
            else:
                sum = sum*math.factorial(i)*math.factorial(n-i)/math.factorial(n)
                coord = (coord - sum) / (2 * (self.level + self.liealg.dual_coxeter()))
            ret_val.append(coord)

        return ret_val

    def get_norm_sym_divisor_ray(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle and normalizes the
        vector by clearing denominators.  
        **DOES NOT WORK WELL WITH FP ARITHMETIC**
        **DOES NOT WORK IN SAGE**

        :return: A list of numbers: the divisor ray given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """
        divisor = self.get_symmetrized_divisor()

        if self.liealg.exact:
            denom_lcm = functools.reduce(lambda x, y: self._lcm(x, y), [long(q.denominator) for q in divisor])
            denom_clear = [long(round(q * denom_lcm)) for q in divisor]
            div_gcd = functools.reduce(lambda x, y: fractions.gcd(x, y), denom_clear)
            if div_gcd > 0:
                return [x//div_gcd for x in denom_clear]
            else:
                return denom_clear
        else:
            n_fact = math.factorial(len(self.weights))
            int_div = [long(round(n_fact * x)) for x in divisor]
            div_gcd = functools.reduce(lambda x, y: fractions.gcd(x, y), int_div)
            if div_gcd > 0:
                return [x // div_gcd for x in int_div]
            else:
                return [x for x in int_div]

    def _lcm(self, x, y):
        return x*y//fractions.gcd(x, y)

    def get_F_curves(self):
        """
        Generates a list of all F-curves with the same number of points as the conformal
        blocks bundle.

        :return: A list of partitions of [1, 2,..., n]: the list of F-curves
        """
        n = len(self.weights)
        all_points = set([x for x in range(1, n+1)])
        ret_list = []
        if n == 3:
            return ret_list

        for r_1 in range(1, n - 2):
            for sset_1 in itertools.combinations(all_points, r_1):
                comp_sset_1 = all_points.difference(sset_1)
                for r_2 in range(1, n - r_1 - 1):
                    for sset_2 in itertools.combinations(comp_sset_1, r_2):
                        comp_sset_2 = comp_sset_1.difference(sset_2)
                        for r_3 in range(1, n - r_1 - r_2):
                            for sset_3 in itertools.combinations(comp_sset_2, r_3):
                                sset_4 = comp_sset_2.difference(sset_3)
                                ret_list.append([sset_1, sset_2, sset_3, tuple(sset_4)])

        return ret_list

    def intersect_F_curve(self, partition):
        """
        Computes the intersection of the divisor associated to this conformal blocks bundle with
        the given F-curve.

        :param partition: A list of 4 lists of integers partitioning the set {1, ..., # points}: the
            F-curve to be intersected.
        :return: An integer: the intersection number.
        """
        ret_val = 0
        wt_list1 = [self.weights[point - 1] for point in partition[0]]
        wt_list2 = [self.weights[point - 1] for point in partition[1]]
        wt_list3 = [self.weights[point - 1] for point in partition[2]]
        wt_list4 = [self.weights[point - 1] for point in partition[3]]

        prod1 = self.liealg.multi_fusion(wt_list1, self.level)
        prod2 = self.liealg.multi_fusion(wt_list2, self.level)
        prod3 = self.liealg.multi_fusion(wt_list3, self.level)
        prod4 = self.liealg.multi_fusion(wt_list4, self.level)

        for wt1 in prod1.keys():
            if prod1[wt1] == 0: continue
            for wt2 in prod2.keys():
                if prod2[wt2] == 0: continue
                for wt3 in prod3.keys():
                    if prod3[wt3] == 0: continue
                    mu_list = [wt1, wt2, wt3]
                    mu_prod = self.liealg.multi_fusion(mu_list, self.level)
                    for wt4 in prod4.keys():
                        if prod4[wt4] == 0: continue
                        if mu_prod[self.liealg.get_dual_weight(wt4)] == 0: continue
                        ret_val += self._degree(wt1, wt2, wt3, wt4, self.level) * \
                                   prod1[wt1] * prod2[wt2] * prod3[wt3] * prod4[wt4]

        return ret_val

    def _degree(self, wt1, wt2, wt3, wt4, level):
        """
        Computes the degree of a four-point conformal blocks vector bundle.  Implements Fakhruddin's
        formula.

        :param wt1: A tuple of integers: a weight of the bundle.
        :param wt2: A tuple of integers: a weight of the bundle.
        :param wt3: A tuple of integers: a weight of the bundle.
        :param wt4: A tuple of integers: a weight of the bundle.
        :param level: A positive integer: the level of the bundle.
        :return: A positive integer: the degree of the bundle.
        """
        liealg = self.liealg
        ret_val = self._get_rank([wt1, wt2, wt3, wt4], level) * (
            liealg.casimirScalar(wt1) + liealg.casimirScalar(wt2) + liealg.casimirScalar(wt3) + liealg.casimirScalar(wt4))

        sum = 0
        prod1 = liealg.fusion(wt1, wt2, level)
        prod2 = liealg.fusion(wt3, wt4, level)
        for mu in prod1.keys():
            mu_star = liealg.get_dual_weight(mu)
            if mu_star in prod2:
                sum += liealg.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        prod1 = liealg.fusion(wt1, wt3, level)
        prod2 = liealg.fusion(wt2, wt4, level)
        for mu in prod1.keys():
            mu_star = liealg.get_dual_weight(mu)
            if mu_star in prod2:
                sum += liealg.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        prod1 = liealg.fusion(wt1, wt4, level)
        prod2 = liealg.fusion(wt2, wt3, level)
        for mu in prod1.keys():
            mu_star = liealg.get_dual_weight(mu)
            if mu_star in prod2:
                sum += liealg.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        ret_val -= sum

        if liealg.exact:
            ret_val = Fraction(ret_val, (2 * (level + liealg.dual_coxeter())))
        else:
            ret_val = round(ret_val / (2 * (level + liealg.dual_coxeter())))

        return ret_val


class SymmetricConformalBlocksBundle(ConformalBlocksBundle):
    """
    A class representing a symmetric conformal blocks vector bundle.
    """

    def __init__(self, liealg, wt, num_points, level):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param wt: A list of integers: the weight of the conformal blocks bundle, repeated at each
            point.
        :param num_points: A positive integer: the number of points of the conformal blocks bundle.
        :param level: A positive integer: the level of the conformal blocks bundle.
        """
        ConformalBlocksBundle.__init__(self, liealg, [wt for i in range(num_points)], level)

    def get_symmetrized_divisor(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle.  Algorithm is
        optimized for the symmetric case.

        :return: A list of numbers: the divisor given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """
        ret_val = []
        n = len(self.weights)
        wt = self.weights[0]
        for i in range(2, n // 2 + 1):
            if self.liealg.exact:
                coord = Fraction(i * (n - i) * self.get_rank() * self.liealg.casimirScalar(wt), n - 1)
            else:
                coord = i * (n - i) * self.get_rank() * self.liealg.casimirScalar(wt) / (n - 1)
            sum_list = [0]
            self._weighted_factor(wt, wt, 1, i - 1, n - i, sum_list, {})

            if self.liealg.exact:
                coord = Fraction(coord - sum_list[0], 2 * (self.level + self.liealg.dual_coxeter()))
            else:
                coord = (coord - sum_list[0]) / (2 * (self.level + self.liealg.dual_coxeter()))
            ret_val.append(coord)

        return ret_val

    def _weighted_factor(self, wt, wt2, mult, wts_rem, ic, ret_val, rank_dict):
        prod = self.liealg.fusion(wt, wt2, self.level)

        for wt3 in prod.keys():
            if wts_rem > 1:
                self._weighted_factor(wt, wt3, mult * prod[wt3], wts_rem - 1, ic, ret_val, rank_dict)
            else:
                if not wt3 in rank_dict:
                    wt_list = [wt for i in range(ic)]
                    wt_list.append(wt3)
                    rank_dict[wt3] = self._get_rank(wt_list, self.level)

                ret_val[0] += self.liealg.casimirScalar(self.liealg.get_dual_weight(wt3)) * mult * prod[wt3] * \
                              rank_dict[wt3]

    def get_sym_F_curves(self):
        """
        Generates a list of all F-curves with the same number of points as the conformal
        blocks bundle, up to permutation of points.

        :return: A list of partitions of [1, 2,..., n]: the list of F-curves
        """
        n = len(self.weights)
        partitions = []

        for part1 in range(int(math.ceil(n / 4)), n - 2):
            for part2 in range(int(math.ceil((n - part1) / 3)), min(n - part1 - 2, part1) + 1):
                for part3 in range(int(math.ceil((n - part1 - part2) / 2)), min(n - part1 - part2 - 1, part2) + 1):
                    part4 = n - part1 - part2 - part3
                    partitions.append((part1, part2, part3, part4))

        ret_list = []
        for partition in partitions:
            p1, p2, p3, p4 = partition[0], partition[1], partition[2], partition[3]
            f_curve = [tuple([x for x in range(1, p1 + 1)]), tuple([x for x in range(p1 + 1, p1 + p2 + 1)]),
                       tuple([x for x in range(p1 + p2 + 1, p1 + p2 + p3 + 1)]),
                       tuple([x for x in range(p1 + p2 + p3 + 1, n + 1)])]
            ret_list.append(f_curve)



        return ret_list