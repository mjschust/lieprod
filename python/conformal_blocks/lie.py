from __future__ import division
from collections import defaultdict
import math, fractions, itertools
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


class SimpleLieAlgebra(object):
    """
    A template class for a simple Lie Algebra.  Objects of this class should be constructed by creating
    an object of the appropriate subclass.  Unimplemented methods must be implemented by
    subclasses of each type.
    """

    def __init__(self, rank, store_fusion=True, exact=True):
        """
        :param rank: A positive integer: the rank of the Lie algebra.
        :param store_fusion: A boolean: if true the lie algebra will save computed fusion products;
        """
        self.rank = rank
        self.store_fusion = store_fusion
        self.exact = exact
        if store_fusion: self._fusion_dict = {}
        self._pos_roots = []
        self._rep_dim_dict = {}
        #self._fte_dict = {}

    def get_rep_dim(self, high_weight):
        """
        Computes the dimension of the representation with given highest weight.
        Implements Weyl's dimension formula.

        :param high_weight: A tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :return: An integer: the dimension of the representation.
        """
        if high_weight in self._rep_dim_dict: return self._rep_dim_dict[high_weight]

        lam = high_weight
        rho = self.get_rho()
        pos_roots = self.get_positive_roots()

        numer = 1
        denom = 1
        for root in pos_roots:
            a = self.killing_form(lam, root)
            b = self.killing_form(rho, root)
            numer = numer * (a + b)
            denom = denom * b


        if self.exact:
            self._rep_dim_dict[high_weight] = Fraction(numer, denom)
        else:
            self._rep_dim_dict[high_weight] = round(numer / denom)
        return self._rep_dim_dict[high_weight]

    def get_dominant_character(self, high_weight):
        """
        Computes the dominant character of the representation.  Implements Freudenthal's recursion
        formula to compute weight multiplicities.  Implementation is heavily influenced by
        the implementation in the LiE system.

        :param high_weight: A tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :return: A dictionary where the keys are tuples and the values are positive integers
            corresponding to the multiplicity of the wt space.
        """

        #Construct the root-level dictionary
        pos_roots = self.get_positive_roots()
        root_level_dict = {}
        for root in pos_roots:
            level = root.get_root_level()
            if level in root_level_dict:
                root_level_dict[level].append(root)
            else:
                root_level_dict[level] = [root]

        #Construct the set of dominant characters
        level = 0
        weight_level_dict = {}
        weight_level_dict[0] = [high_weight]
        dom_weights = set()
        dom_weights.add(high_weight)
        while True:
            done = True
            for key in weight_level_dict.keys():
                if level <= key:
                    done = False
                    break
            if done:
                break

            if not level in weight_level_dict:
                level = level + 1
                continue

            for wt in weight_level_dict[level]:
                for root_lev in root_level_dict.keys():
                    for root in root_level_dict[root_lev]:
                        new_weight = self._sub_weights(wt, root)
                        if self.is_dominant(new_weight):
                            if level + root_lev in weight_level_dict:
                                if not new_weight in weight_level_dict[level + root_lev]:
                                    weight_level_dict[level + root_lev].add(new_weight)
                                    dom_weights.add(new_weight)
                            else:
                                weight_level_dict[level + root_lev] = {new_weight}
                                dom_weights.add(new_weight)

            level = level + 1

        #Calculate multiplicities of dominant weights
        dom_char = {}
        sorted_levels = sorted(weight_level_dict.keys())
        for level in sorted_levels:
            for wt in weight_level_dict[level]:
                self._compute_mult(high_weight, wt, pos_roots, dom_weights, dom_char)

        return dom_char

    def _compute_mult(self, high_weight, wt, pos_roots, dom_weights, dom_char):
        """
        This implements Freudenthal's recursion formula.  Expects a dominant wt in _dom_weights.
        """
        if wt in dom_char:
            return dom_char[wt]
        if wt == high_weight:
            dom_char[wt] = 1
            return 1

        mult_sum = 0
        for root in pos_roots:
            n = 0
            new_weight = wt
            a = self.killing_form(wt, root)
            b = self.killing_form(root, root)
            while True:
                n = n + 1
                new_weight = self._add_weights(new_weight, root)
                new_dom_weight = self.reflect_to_chamber(new_weight)
                if not new_dom_weight in dom_weights: break

                mult_sum = mult_sum + (a + n * b) * self._compute_mult(high_weight, new_dom_weight, pos_roots, dom_weights, dom_char)

        rho = self.get_rho()
        if self.exact:
            multiplicity = Fraction(2 * mult_sum, (self.length_squared(self._add_weights(high_weight, rho)) - self.length_squared(
                    self._add_weights(wt, rho))))
        else:
            multiplicity = 2 * mult_sum / (
                    self.length_squared(self._add_weights(high_weight, rho)) - self.length_squared(self._add_weights(wt, rho)))
        dom_char[wt] = multiplicity
        return multiplicity

    def tensor(self, wt1, wt2):
        """
        Computes the tensor product decomposition of the irreducible representations with highest
        weights wt1 and wt2.

        :param wt1: A tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :return: A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the tensor product of wt1 and wt2.
        """

        # Want wt1 to have larger dimension
        if self.get_rep_dim(wt1) < self.get_rep_dim(wt2):
            wt1, wt2 = wt2, wt1

        rho = self.get_rho()
        dom_char = self.get_dominant_character(wt2)
        lam_rho_sum = self._add_weights(wt1, rho)
        ret_dict = {}

        # Traverse entire character
        for dom_weight in dom_char.keys():
            for orbit_weight in self.get_orbit_iter(dom_weight):
                new_sum = self._add_weights(lam_rho_sum, orbit_weight)
                new_dom_weight, parity = self.reflect_to_chamber_with_parity(new_sum)
                new_dom_weight = self._sub_weights(new_dom_weight, rho)
                if not self.is_dominant(new_dom_weight): continue
                if new_dom_weight in ret_dict:
                    ret_dict[new_dom_weight] = ret_dict[new_dom_weight] + dom_char[dom_weight] * parity
                else:
                    ret_dict[new_dom_weight] = dom_char[dom_weight] * parity

        return ret_dict

    def fusion(self, wt1, wt2, ell):
        """
        Computes the fusion product decomposition of the irreducible representations with highest
        weights wt1 and wt2 at level ell.

        :param wt1: A tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :param ell: A positive integer: corresponds to the level of the fusion product.
        :return:  A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the fusion product of wt1 and wt2.
        """

        if self.store_fusion and (wt1, wt2, ell) in self._fusion_dict:
            return self._fusion_dict[(wt1, wt2, ell)]

        ten_decom = self.tensor(wt1, wt2)
        ret_dict = {}
        rho = self.get_rho()

        for wt in ten_decom.keys():
            if self.get_level(wt) == ell + 1: continue

            wt_rho = self._add_weights(wt, rho)
            new_weight, parity = self.reflect_to_alcove_with_parity(wt_rho, ell + self.get_level(rho) + 1)
            lev_ell_weight = self._sub_weights(new_weight, rho)
            if not self.is_dominant(lev_ell_weight) or self.get_level(lev_ell_weight) > ell: continue

            if lev_ell_weight in ret_dict:
                ret_dict[lev_ell_weight] = ret_dict[lev_ell_weight] + ten_decom[wt] * parity
            else:
                ret_dict[lev_ell_weight] = ten_decom[wt] * parity

        if self.store_fusion:
            self._fusion_dict[(wt1, wt2, ell)] = ret_dict

        return ret_dict

    def multi_fusion(self, wts, level):
        """
        Computes the fusion product of a list of representations.

        :param wts: A list of tuples of integers: the list of weights.
        :param level: A positive integer: corresponds to the level of the fusion product.
        :return: A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the fusion product of the weights.
        """
        wt_dict = defaultdict(int)
        rem_wts = list(wts)
        cur_wt = rem_wts.pop()
        wt_dict[cur_wt] = 1
        while (len(rem_wts) > 0):
            new_wt_dict = defaultdict(int)
            cur_wt = rem_wts.pop()
            for wt in wt_dict.keys():
                prod = self.fusion(cur_wt, wt, level)
                for wt2 in prod.keys():
                    new_wt_dict[wt2] += wt_dict[wt] * prod[wt2]

            wt_dict = new_wt_dict

        return wt_dict


    def killing_form(self, wt1, wt2):
        """
        Computes the Killing form product of two weights.

        :param wt1: A tuple of numbers: a weight.
        :param wt2: A tuple of numbers: a weight.
        :return: A number: the Killing form product of wt1 and wt2.
        """
        raise NotImplementedError

    def length_squared(self, wt):
        """
        Computes the squared length of a weight using the Killing form of the algebra.

        :param wt: A tuple of numbers: a weight.
        :return: A positive number: the squared length of wt.
        """
        return self.killing_form(wt, wt)

    def casimirScalar(self, wt):
        """
        Computes the Casimir scalar of a weight.

        :param wt: A tuple of numbers: a weight.
        :return: A number: the Casimir scalar of wt.
        """
        twoRho = tuple([2 for i in range(self.rank)])
        wt2 = self._add_weights(wt, twoRho)
        return self.killing_form(wt, wt2)

    def dual_coxeter(self):
        """
        Computes the dual coxeter number of the Lie algebra.

        :return: An integer: the dual coxeter number.
        """
        raise NotImplementedError

    def is_dominant(self, wt):
        """
        Checks if the weight is dominant.

        :param wt: A tuple of numbers; a weight of the Lie algebra.
        :return: A boolean.
        """
        for coord in wt:
            if coord < 0:
                return False

        return True

    def get_level(self, wt):
        """
        Computes the level of the weight, which is defined as the Killing product of wt and
        the highest root of the algebra.

        :param wt:  A tuple of numbers: a weight.
        :return: A number: the level of the weight.
        """
        raise NotImplementedError

    def get_dual_weight(self, wt):
        """
        Computes the highest weight of the contragredient representation associated to wt.

        :param wt: A tuple of positive integers: a dominant integral weight.
        :return: A tuple of positive integers: a dominant integral weight.
        """
        raise NotImplementedError

    def get_rho(self):
        """
        Computes one-half the sum of the positive roots of the algebra, which is usually denoted
        by rho.

        :return: A tuple of positive integers: the weight rho.
        """
        ret_coords = []

        for i in range(self.rank):
            ret_coords.append(1)

        return tuple(ret_coords)

    def get_positive_roots(self):
        """
        Computes a list of all positive roots of the Lie algebra.

        :return: A tuple of weights: the positive roots of the algebra.
        """
        raise NotImplementedError

    def get_weights(self, level):
        """
        Computes a list of all weights of level less than or equal to the given level.

        :param level: A positive integer.
        :return: A list of weights with level less than level.
        """
        raise NotImplementedError

    def reflect_to_chamber(self, wt):
        """
        Reflects a weight into the dominant chamber of the the Lie algebra.

        :param wt: A tuple of numbers: a weight.
        :return: A tuple of non-negative numbers: a dominant weight.
        """
        raise NotImplementedError

    def reflect_to_chamber_with_parity(self, wt):
        """
        Reflects a weight into the dominant chamber of the the Lie algebra, keeping track of the
        parity of the reflections.

        :param wt: A tuple of numbers: a weight.
        :return: A tuple of non-negative numbers: a dominant weight; and a number equal to +1 or -1.
        """
        raise NotImplementedError

    def reflect_to_alcove_with_parity(self, wt, ell):
        """
        Reflects a weight into the level ell fundamental chamber of the the Lie algebra, keeping track of the
        parity of the reflections.

        :param wt: A tuple of numbers: a weight.
        :param ell: A positive integer: the level.
        :return: A tuple of non-negative numbers: a dominant weight; and a number equal to +1 or -1.
        """
        raise NotImplementedError

    def get_orbit_iter(self, wt):
        """
        Returns an iterable object that iterates through the Weyl group orbit of the given weight.

        :param wt: A tuple of numbers: a weight.
        :return: An iterator object.
        """
        raise NotImplementedError

    def _convert_funds_to_epsilons(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_epsilons_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_funds_to_roots(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_roots_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def _add_weights(self, wt1, wt2):
        '''
        Adds two weights and returns the sum as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] + wt2[i])

        return tuple(ret_coords)

    def _sub_weights(self, wt1, wt2):
        '''
        Subtracts wt1 from wt2 the difference as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] - wt2[i])

        return tuple(ret_coords)



class TypeALieAlgebra(SimpleLieAlgebra):
    """
    A type A Lie algebra.
    """

    def killing_form(self, wt1, wt2):
        ret_val = 0
        ep_coords1 = self._convert_funds_to_epsilons(wt1)
        ep_coords2 = self._convert_funds_to_epsilons(wt2)

        sum1 = sum2 = 0
        for i in range(self.rank + 1):
            ret_val += ep_coords1[i] * ep_coords2[i]
            sum1 += ep_coords1[i]
            sum2 += ep_coords2[i]

        if self.exact:
            ret_val -= Fraction(sum1 * sum2, self.rank + 1)
        else:
            ret_val -= sum1 * sum2 / (self.rank + 1)

        return ret_val

    def dual_coxeter(self):
        return self.rank + 1

    def get_level(self, wt):
        return sum(wt)

    def get_dual_weight(self, wt):
        return tuple(wt[::-1])

    def get_positive_roots(self):
        if len(self._pos_roots) > 0: return self._pos_roots

        ret_list = []
        coords = []
        for i in range(self.rank):
            coords.append(0)

        for i in range(len(coords)):
            for j in range(i, len(coords)):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        self._pos_roots = ret_list
        return ret_list

    def get_weights(self, level):
        """
        Computes a list of all weights of level less than or equal to the given level.

        :param level: A positive integer.
        :return: A list of weights with level less than level.
        """
        return [tuple(coords) for coords in self._get_weights(level, self.rank)]

    def _get_weights(self, level, rank):
        ret_list = []
        if rank == 1:
            for i in range(level + 1):
                ret_list.append([i])
        else:
            r_minus_one_list = self._get_weights(level, rank - 1)
            for coord in r_minus_one_list:
                for i in range(level - sum(coord) + 1):
                    ret_list.append(coord + [i])

        return ret_list

    def reflect_to_chamber(self, wt):
        ret_coords = self._insertsort(self._convert_funds_to_epsilons(wt))

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return self._convert_epsilons_to_funds(ret_coords)

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return self._convert_epsilons_to_funds(ret_coords), parity

    def _insertsort_parity(self, coords):
        ret_list = list(coords)

        parity = 1
        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                parity = parity * -1
                j = j - 1

        return ret_list, parity

    def reflect_to_alcove_with_parity(self, wt, ell):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))
        ret_coords = [x - ret_coords[-1] for x in ret_coords]

        while (ret_coords[0] > ell):
            ret_coords[-1] = ret_coords[0] - ell
            ret_coords[0] = ell
            ret_coords, fin_parity = self._insertsort_parity(ret_coords)
            ret_coords = [x - ret_coords[-1] for x in ret_coords]
            parity = parity * -1 * fin_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def get_orbit_iter(self, wt):
        return self._TypeAOrbitIterator(self, wt)

    def _convert_funds_to_epsilons(self, coords):
        #if coords in self._fte_dict: return self._fte_dict[coords]

        ret_coords = [0]
        part = 0
        for i in reversed(range(len(coords))):
            part += coords[i]
            #Looks bad but faster than using a deque for usual case rank
            ret_coords.insert(0, part)

        #self._fte_dict[coords] = ret_coords
        return ret_coords

    def _convert_epsilons_to_funds(self, coords):

        ret_coords = []
        for i in range(len(coords) - 1):
            ret_coords.append(coords[i] - coords[i + 1])

        return tuple(ret_coords)

    def _convert_roots_to_funds(self, coords):
        if len(coords) == 1: return [2 * coords[0]]

        ret_coords = []
        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 1):
            ret_coords.append(2 * coords[i] - coords[i + 1] - coords[i - 1])

        ret_coords.append(2 * coords[-1] - coords[-2])
        return ret_coords

    class _TypeAOrbitIterator(object):
        '''
        Optimized iterator object that traverses the Weyl group orbit of a given weight.

        Attributes:
            no public attributes
        '''
        def __init__(self, liealg, wt):
            ep_coords = liealg._convert_funds_to_epsilons(liealg.reflect_to_chamber(wt))

            #Construct list of items and multiplicities
            self._item_list = [ep_coords[0]]
            cur_item = ep_coords[0]
            rem_list = [0]
            for item in ep_coords:
                if item < cur_item:
                    self._item_list.append(item)
                    rem_list.append(1)
                    cur_item = item
                else:
                    rem_list[-1] += 1

            #Contruct matrix of remaining items, and initial index list
            index_list = []
            rem_mat = [list(rem_list)]
            for i in range(len(ep_coords)):
                j = 0
                while rem_mat[i][j] == 0:
                    j += 1
                index_list.append(j)
                rem_mat.append(list(rem_mat[i]))
                rem_mat[i+1][j] -= 1

            self._index_list = index_list
            self._rem_mat = rem_mat
            self.done = False
            self.liealg = liealg

        def __iter__(self):
            return self

        def next(self):
            if self.done: raise StopIteration()
            r = len(self._index_list)
            num_items = len(self._item_list)

            #Construct new weight
            ep_coords = []
            for index in self._index_list:
                ep_coords.append(self._item_list[index])
            ret_val = self.liealg._convert_epsilons_to_funds(ep_coords)

            #Find index to increment
            i = r-2
            j = 0
            while i >= 0:
                j = self._index_list[i] + 1
                while j < num_items:
                    if self._rem_mat[i][j] > 0: break
                    j += 1
                if j < num_items: break
                i -= 1

            #If we're finished, return the last weight
            if i < 0:
                self.done = True
                return ret_val

            #Increment indices
            self._index_list[i] = j
            self._rem_mat[i+1] = list(self._rem_mat[i])
            self._rem_mat[i+1][j] -= 1
            i += 1
            while i < r:
                j = 0
                while self._rem_mat[i][j] == 0:
                    j += 1
                self._index_list[i] = j
                self._rem_mat[i + 1] = list(self._rem_mat[i])
                self._rem_mat[i + 1][j] -= 1
                i += 1

            return ret_val

class TypeBLieAlgebra(SimpleLieAlgebra):
    """
    A type B Lie algebra.
    """

    def __init__(self, rank, **kwargs):
        """
        :param rank: A positive integer: the rank of the Lie algebra.
        :param store_fusion: A boolean: if true the lie algebra will save computed fusion products;
        """
        if rank < 2:
            raise ValueError("Lie Algebra does not exist.")

        super(TypeBLieAlgebra, self).__init__(rank, **kwargs)

    def killing_form(self, wt1, wt2):
        ret_val = 0
        ep_coords1 = self._convert_funds_to_epsilons(wt1)
        ep_coords2 = self._convert_funds_to_epsilons(wt2)

        for i in range(self.rank):
            ret_val = ret_val + ep_coords1[i] * ep_coords2[i]

        return ret_val

    def dual_coxeter(self):
        return 2*self.rank - 1

    def get_level(self, wt):
        if len(wt) == 2:
            return wt[0] + wt[1]
        else:
            ret_val = wt[0] + wt[-1]
            for i in range(1, len(wt)-1):
                ret_val += 2*wt[i]
            return ret_val

    def get_dual_weight(self, wt):
        return tuple(wt)

    def get_positive_roots(self):
        if len(self._pos_roots) > 0: return self._pos_roots

        ret_list = []
        coords = []
        for i in range(self.rank):
            coords.append(0)

        for i in range(len(coords)):
            for j in range(i, len(coords)):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        for i in range(len(coords)-1, 0, -1):
            coords[i] = 2
            for j in range(i-1, -1, -1):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i-1, -1, -1):
                coords[j] = 0

        self._pos_roots = ret_list
        return ret_list

    def get_weights(self, level):
        ret_list = []
        for a_1 in range(level + 1):
            coords_list = self._get_weights(level - a_1, self.rank - 1)
            for coords in coords_list:
                ret_list.append(tuple(itertools.chain([a_1], coords)))
        return ret_list

    def _get_weights(self, level, rank):
        ret_list = []
        if rank == 1:
            for i in range(level + 1):
                ret_list.append([i])
        else:
            for a_i in range(level//2 + 1):
                coords_list = self._get_weights(level - 2*a_i, rank - 1)
                for coords in coords_list:
                    ret_list.append([a_i] + coords)

        return ret_list

    def reflect_to_chamber(self, wt):
        #First reflect by making epsilon coords positive
        ret_coords = self._convert_funds_to_epsilons(wt)

        ret_coords = [abs(x) for x in ret_coords]

        #Sort to finish reflection into chamber
        ret_coords = self._insertsort(ret_coords)

        #Return weight in terms of fundamental weights
        return self._convert_epsilons_to_funds(ret_coords)

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        #First reflect by making epsilon coords positive
        ret_coords = self._convert_funds_to_epsilons(wt)
        parity = 1
        for i in range(len(ret_coords)):
            if ret_coords[i] < 0:
                ret_coords[i] = -ret_coords[i]
                parity *= -1

        #Then sort to finish reflection
        ret_coords, sort_parity = self._insertsort_parity(ret_coords)
        parity *= sort_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def _insertsort_parity(self, coords):
        ret_list = list(coords)

        parity = 1
        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                parity = parity * -1
                j = j - 1

        return ret_list, parity

    def reflect_to_alcove_with_parity(self, wt, ell):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))

        while (ret_coords[0] + ret_coords[1] > ell):
            #wt := wt + (ell-level(wt))*theta
            ret_coords[0], ret_coords[1] = ell - ret_coords[1], ell - ret_coords[0]

            #Return to chamber
            fin_parity = -1
            for i in range(len(ret_coords)):
                if ret_coords[i] < 0:
                    ret_coords[i] = -ret_coords[i]
                    fin_parity *= -1

            ret_coords, sort_parity = self._insertsort_parity(ret_coords)
            fin_parity *= sort_parity

            parity *= fin_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def get_orbit_iter(self, wt):
        return self._TypeBOrbitIterator(self, wt)

    def _convert_funds_to_epsilons(self, coords):
        #if coords in self._fte_dict: return list(self._fte_dict[coords])

        ret_coords = []
        if self.exact:
            part = Fraction(coords[-1], 2)
        else:
            part = coords[-1] / 2
        ret_coords.append(part)
        for i in reversed(range(len(coords)-1)):
            part += coords[i]
            # Looks bad but faster than using a deque for usual case rank
            ret_coords.insert(0, part)

        #self._fte_dict[coords] = ret_coords
        return ret_coords

    def _convert_epsilons_to_funds(self, coords):
        ret_coords = []
        for i in range(len(coords) - 1):
            ret_coords.append(coords[i] - coords[i + 1])
        ret_coords.append(2 * coords[-1])

        return tuple(ret_coords)

    def _convert_roots_to_funds(self, coords):
        if len(coords) == 2:
            return [2 * coords[0] - coords[1], -2 * coords[0] + 2 * coords[1]]

        ret_coords = []
        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 1):
            ret_coords.append(2 * coords[i] - coords[i + 1] - coords[i - 1])

        ret_coords.append(2 * coords[-1] - 2 * coords[-2])
        return tuple(ret_coords)

    class _TypeBOrbitIterator(object):
        '''
        Optimized iterator object that traverses the Weyl group orbit of a given weight.

        Attributes:
            no public attributes
        '''
        def __init__(self, liealg, wt):
            ep_coords = liealg._convert_funds_to_epsilons(liealg.reflect_to_chamber(wt))

            #Construct list of items and multiplicities
            self._item_list = [ep_coords[0]]
            cur_item = ep_coords[0]
            rem_list = [0]
            for item in ep_coords:
                if item < cur_item:
                    self._item_list.append(item)
                    rem_list.append(1)
                    cur_item = item
                else:
                    rem_list[-1] += 1

            #Contruct matrix of remaining items, and initial index list
            index_list = []
            rem_mat = [list(rem_list)]
            for i in range(len(ep_coords)):
                j = 0
                while rem_mat[i][j] == 0:
                    j += 1
                index_list.append(j)
                rem_mat.append(list(rem_mat[i]))
                rem_mat[i+1][j] -= 1

            self._index_list = index_list
            self._rem_mat = rem_mat
            self.perms_done = False
            self.liealg = liealg

            # Get first permutation, and construct index iterator, which iterates over subsets of
            # the indices of the non-zero elements of the permutation
            self.cur_perm = self._next_perm()
            non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
            self._index_iter = itertools.chain.from_iterable(itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds) + 1))

        def __iter__(self):
            return self

        def next(self):
            try:
                neg_inds = self._index_iter.next()
                ep_coords = list(self.cur_perm)
            except StopIteration:
                if self.perms_done:
                    raise StopIteration()
                else:
                    self.cur_perm = self._next_perm()
                    ep_coords = list(self.cur_perm)
                    non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
                    self._index_iter = itertools.chain.from_iterable(itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds)+1))
                    neg_inds = self._index_iter.next()

            for i in neg_inds:
                ep_coords[i] = -ep_coords[i]

            return self.liealg._convert_epsilons_to_funds(ep_coords)

        def _next_perm(self):
            r = len(self._index_list)
            num_items = len(self._item_list)

            # Construct new weight
            ep_coords = []
            for index in self._index_list:
                ep_coords.append(self._item_list[index])

            # Find index to increment
            i = r - 2
            j = 0
            while i >= 0:
                j = self._index_list[i] + 1
                while j < num_items:
                    if self._rem_mat[i][j] > 0: break
                    j += 1
                if j < num_items: break
                i -= 1

            # If we're finished, return the last weight
            if i < 0:
                self.perms_done = True
                return ep_coords

            # Increment indices
            self._index_list[i] = j
            self._rem_mat[i + 1] = list(self._rem_mat[i])
            self._rem_mat[i + 1][j] -= 1
            i += 1
            while i < r:
                j = 0
                while self._rem_mat[i][j] == 0:
                    j += 1
                self._index_list[i] = j
                self._rem_mat[i + 1] = list(self._rem_mat[i])
                self._rem_mat[i + 1][j] -= 1
                i += 1

            return ep_coords


class TypeCLieAlgebra(SimpleLieAlgebra):
    """
    A type C Lie algebra.
    """

    def __init__(self, rank, **kwargs):
        """
        :param rank: A positive integer: the rank of the Lie algebra.
        :param store_fusion: A boolean: if true the lie algebra will save computed fusion products;
        """
        if rank < 2:
            raise ValueError("Lie Algebra does not exist.")

        super(TypeCLieAlgebra, self).__init__(rank, **kwargs)

    def killing_form(self, wt1, wt2):
        ret_val = 0
        ep_coords1 = self._convert_funds_to_epsilons(wt1)
        ep_coords2 = self._convert_funds_to_epsilons(wt2)

        for i in range(self.rank):
            ret_val = ret_val + ep_coords1[i] * ep_coords2[i]

        if self.exact:
            return Fraction(ret_val, 2)
        else:
            return ret_val/2

    def dual_coxeter(self):
        return self.rank + 1

    def get_level(self, wt):
        return sum(wt)

    def get_dual_weight(self, wt):
        return tuple(wt)

    def get_positive_roots(self):
        if len(self._pos_roots) > 0: return self._pos_roots

        ret_list = []
        coords = []
        for i in range(self.rank):
            coords.append(0)

        for i in range(len(coords)):
            for j in range(i, len(coords)):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        coords[-1] = 1
        for i in range(len(coords)-2, -1, -1):
            coords[i] = 2
            ret_list.append(_Root(self, list(coords)))
            for j in range(i-1, -1, -1):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i-1, -1, -1):
                coords[j] = 0

        self._pos_roots = ret_list
        return ret_list

    def get_weights(self, level):
        """
        Computes a list of all weights of level less than or equal to the given level.

        :param level: A positive integer.
        :return: A list of weights with level less than level.
        """
        return [tuple(coords) for coords in self._get_weights(level, self.rank)]

    def _get_weights(self, level, rank):
        ret_list = []

        r_minus_one_list = self._get_weights(level, rank - 1)
        for coord in r_minus_one_list:
            for i in range(level - sum(coord) + 1):
                ret_list.append(coord + [i])

        return ret_list

    def reflect_to_chamber(self, wt):
        #First reflect by making epsilon coords positive
        ret_coords = self._convert_funds_to_epsilons(wt)

        ret_coords = [abs(x) for x in ret_coords]

        #Sort to finish reflection into chamber
        ret_coords = self._insertsort(ret_coords)

        #Return weight in terms of fundamental weights
        return self._convert_epsilons_to_funds(ret_coords)

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        #First reflect by making epsilon coords positive
        ret_coords = self._convert_funds_to_epsilons(wt)
        parity = 1
        for i in range(len(ret_coords)):
            if ret_coords[i] < 0:
                ret_coords[i] = -ret_coords[i]
                parity *= -1

        #Then sort to finish reflection
        ret_coords, sort_parity = self._insertsort_parity(ret_coords)
        parity *= sort_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def _insertsort_parity(self, coords):
        ret_list = list(coords)

        parity = 1
        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                parity = parity * -1
                j = j - 1

        return ret_list, parity

    def reflect_to_alcove_with_parity(self, wt, ell):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))

        while (ret_coords[0] > ell):
            #wt := wt + (ell-level(wt))*theta
            ret_coords[0] = 2*ell-ret_coords[0]

            #Return to chamber
            fin_parity = -1
            for i in range(len(ret_coords)):
                if ret_coords[i] < 0:
                    ret_coords[i] = -ret_coords[i]
                    fin_parity *= -1

            ret_coords, sort_parity = self._insertsort_parity(ret_coords)
            fin_parity *= sort_parity

            parity *= fin_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def get_orbit_iter(self, wt):
        return self._TypeCOrbitIterator(self, wt)

    def _convert_funds_to_epsilons(self, coords):
        #if coords in self._fte_dict: return list(self._fte_dict[coords])

        ret_coords = []
        part = 0
        for i in reversed(range(len(coords))):
            part += coords[i]
            # Looks bad but faster than using a deque for usual case rank
            ret_coords.insert(0, part)

        #self._fte_dict[coords] = ret_coords
        return ret_coords

    def _convert_epsilons_to_funds(self, coords):
        ret_coords = []
        for i in range(len(coords) - 1):
            ret_coords.append(coords[i] - coords[i + 1])
        ret_coords.append(coords[-1])

        return tuple(ret_coords)

    def _convert_roots_to_funds(self, coords):
        if len(coords) == 2:
            return [2 * coords[0] - 2 * coords[1], -coords[0] + 2 * coords[1]]

        ret_coords = []
        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 2):
            ret_coords.append(2 * coords[i] - coords[i + 1] - coords[i - 1])

        ret_coords.append(2 * coords[-2] - 2 * coords[-1] - coords[-3])
        ret_coords.append(2 * coords[-1] - coords[-2])
        return tuple(ret_coords)

    class _TypeCOrbitIterator(object):
        '''
        Optimized iterator object that traverses the Weyl group orbit of a given weight.

        Attributes:
            no public attributes
        '''
        def __init__(self, liealg, wt):
            ep_coords = liealg._convert_funds_to_epsilons(liealg.reflect_to_chamber(wt))

            #Construct list of items and multiplicities
            self._item_list = [ep_coords[0]]
            cur_item = ep_coords[0]
            rem_list = [0]
            for item in ep_coords:
                if item < cur_item:
                    self._item_list.append(item)
                    rem_list.append(1)
                    cur_item = item
                else:
                    rem_list[-1] += 1

            #Contruct matrix of remaining items, and initial index list
            index_list = []
            rem_mat = [list(rem_list)]
            for i in range(len(ep_coords)):
                j = 0
                while rem_mat[i][j] == 0:
                    j += 1
                index_list.append(j)
                rem_mat.append(list(rem_mat[i]))
                rem_mat[i+1][j] -= 1

            self._index_list = index_list
            self._rem_mat = rem_mat
            self.perms_done = False
            self.liealg = liealg

            # Get first permutation, and construct index iterator, which iterates over subsets of
            # the indices of the non-zero elements of the permutation
            self.cur_perm = self._next_perm()
            non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
            self._index_iter = itertools.chain.from_iterable(itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds) + 1))

        def __iter__(self):
            return self

        def next(self):
            try:
                neg_inds = self._index_iter.next()
                ep_coords = list(self.cur_perm)
            except StopIteration:
                if self.perms_done:
                    raise StopIteration()
                else:
                    self.cur_perm = self._next_perm()
                    ep_coords = list(self.cur_perm)
                    non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
                    self._index_iter = itertools.chain.from_iterable(itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds)+1))
                    neg_inds = self._index_iter.next()

            for i in neg_inds:
                ep_coords[i] = -ep_coords[i]

            return self.liealg._convert_epsilons_to_funds(ep_coords)

        def _next_perm(self):
            r = len(self._index_list)
            num_items = len(self._item_list)

            # Construct new weight
            ep_coords = []
            for index in self._index_list:
                ep_coords.append(self._item_list[index])

            # Find index to increment
            i = r - 2
            j = 0
            while i >= 0:
                j = self._index_list[i] + 1
                while j < num_items:
                    if self._rem_mat[i][j] > 0: break
                    j += 1
                if j < num_items: break
                i -= 1

            # If we're finished, return the last weight
            if i < 0:
                self.perms_done = True
                return ep_coords

            # Increment indices
            self._index_list[i] = j
            self._rem_mat[i + 1] = list(self._rem_mat[i])
            self._rem_mat[i + 1][j] -= 1
            i += 1
            while i < r:
                j = 0
                while self._rem_mat[i][j] == 0:
                    j += 1
                self._index_list[i] = j
                self._rem_mat[i + 1] = list(self._rem_mat[i])
                self._rem_mat[i + 1][j] -= 1
                i += 1

            return ep_coords

class TypeDLieAlgebra(SimpleLieAlgebra):
    """
    A type D Lie algebra.
    """

    def __init__(self, rank, **kwargs):
        """
        :param rank: A positive integer: the rank of the Lie algebra.
        :param store_fusion: A boolean: if true the lie algebra will save computed fusion products;
        """
        if rank < 3:
            raise ValueError("Lie Algebra does not exist.")

        super(TypeDLieAlgebra, self).__init__(rank, **kwargs)

    def killing_form(self, wt1, wt2):
        ret_val = 0
        ep_coords1 = self._convert_funds_to_epsilons(wt1)
        ep_coords2 = self._convert_funds_to_epsilons(wt2)

        for i in range(self.rank):
            ret_val = ret_val + ep_coords1[i] * ep_coords2[i]

        return ret_val

    def dual_coxeter(self):
        return 2 * self.rank - 2

    def get_level(self, wt):
        ret_val = wt[0]

        for i in range(1, len(wt) - 2):
            ret_val += 2 * wt[i]

        ret_val += wt[-2] + wt[-1]
        return ret_val

    def get_dual_weight(self, wt):
        if self.rank % 2 == 0:
            return tuple(wt)
        else:
            ret_coords = list(wt)
            ret_coords[-2], ret_coords[-1] = ret_coords[-1], ret_coords[-2]
            return tuple(ret_coords)

    def get_positive_roots(self):
        if len(self._pos_roots) > 0:
            return self._pos_roots

        ret_list = []
        coords = [0 for i in range(self.rank)]

        for i in range(self.rank - 2):
            for j in range(i, self.rank - 2):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            coords[-2] = 1
            ret_list.append(_Root(self, list(coords)))
            coords[-2] = 0
            coords[-1] = 1
            ret_list.append(_Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        coords[-2] = 1
        ret_list.append(_Root(self, list(coords)))
        coords[-2] = 0
        coords[-1] = 1
        ret_list.append(_Root(self, list(coords)))
        coords[-2] = 1
        for i in range(len(coords)-3, -1, -1):
            for j in range(i, -1, -1):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i, -1, -1):
                coords[j] = 0

            coords[i] = 2

        self._pos_roots = ret_list
        return ret_list

    def get_weights(self, level):
        """
        Computes a list of all weights of level less than or equal to the given level.

        :param level: A positive integer.
        :return: A list of weights with level less than level.
        """
        ret_list = []
        for i in range(level + 1):
            for j in range(level - i + 1):
                ret_list.extend([tuple(coords + [i, j]) for coords in self._get_weights(level - i - j, self.rank - 2)])
        return ret_list

    def _get_weights(self, level, rank):
        ret_list = []
        if rank == 1:
            for i in range(level + 1):
                ret_list.append([i])
        else:
            r_minus_one_list = self._get_weights(level, rank - 1)
            for coord in r_minus_one_list:
                for i in range((level - coord[0] - 2*sum(coord[1:])) // 2 + 1):
                    ret_list.append(coord + [i])

        return ret_list

    def reflect_to_chamber(self, wt):
        #First reflect by making epsilon coords positive, keeping track of cum sign
        ret_coords = self._convert_funds_to_epsilons(wt)

        sign = 1
        for i in range(len(ret_coords)):
            if ret_coords[i] < 0:
                sign *= -1
            ret_coords[i] = abs(ret_coords[i])

        #Sort coords
        ret_coords = self._insertsort(ret_coords)

        #Multiply smallest coord by cum sign
        ret_coords[-1] *= sign

        #Return weight in terms of fundamental weights
        return self._convert_epsilons_to_funds(ret_coords)

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        ret_coords = self._convert_funds_to_epsilons(wt)
        ret_coords, parity = self._reflect_eps_coords_to_chamber_with_parity(ret_coords)

        return self._convert_epsilons_to_funds(ret_coords), parity

    def _reflect_eps_coords_to_chamber_with_parity(self, coords):
        #First 'reflect' by making epsilon coords positive, keeping track of negative
        #signs
        ret_coords = list(coords)
        sign = 1
        for i in range(len(ret_coords)):
            if ret_coords[i] < 0:
                ret_coords[i] = -ret_coords[i]
                sign *= -1

        #Then sort and multiply last coord by above sign to finish reflection
        ret_coords, parity = self._insertsort_parity(ret_coords)
        ret_coords[-1] *= sign

        return ret_coords, parity

    def _insertsort_parity(self, coords):
        ret_list = list(coords)

        parity = 1
        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                parity = parity * -1
                j = j - 1

        return ret_list, parity

    def reflect_to_alcove_with_parity(self, wt, ell):
        ret_coords = self._convert_funds_to_epsilons(wt)
        ret_coords, parity = self._reflect_eps_coords_to_chamber_with_parity(ret_coords)

        while (ret_coords[0] + ret_coords[1] > ell):
            #wt := wt + (ell-level(wt))*theta
            ret_coords[0], ret_coords[1] = ell - ret_coords[1], ell - ret_coords[0]

            #Return to chamber
            ret_coords, chamber_parity = self._reflect_eps_coords_to_chamber_with_parity(ret_coords)

            parity *= -1 * chamber_parity

        return self._convert_epsilons_to_funds(ret_coords), parity

    def get_orbit_iter(self, wt):
        return self._TypeDOrbitIterator(self, wt)

    def _convert_funds_to_epsilons(self, coords):
        #if coords in self._fte_dict: return list(self._fte_dict[coords])

        ret_coords = []
        if self.exact:
            part = Fraction(-coords[-2] + coords[-1], 2)
            ret_coords.insert(0, part)
            part = Fraction(coords[-2] + coords[-1], 2)
            ret_coords.insert(0, part)
        else:
            part = (-coords[-2] + coords[-1]) / 2
            ret_coords.insert(0, part)
            part = (coords[-2] + coords[-1]) / 2
            ret_coords.insert(0, part)
        for i in reversed(range(len(coords) - 2)):
            part += coords[i]
            # Looks bad but faster than using a deque for usual case rank
            ret_coords.insert(0, part)

        #self._fte_dict[coords] = ret_coords
        return ret_coords

    def _convert_epsilons_to_funds(self, coords):
        ret_coords = []
        for i in range(len(coords) - 2):
            ret_coords.append(coords[i] - coords[i + 1])
        ret_coords.append(coords[-2] - coords[-1])
        ret_coords.append(coords[-2] + coords[-1])

        return tuple(ret_coords)

    def _convert_roots_to_funds(self, coords):
        ret_coords = []

        if self.rank == 3:
            ret_coords.append(2 * coords[-3] - coords[-2] - coords[-1])
            ret_coords.append(-coords[-3] + 2 * coords[-2])
            ret_coords.append(-coords[-3] + 2 * coords[-1])
            return tuple(ret_coords)

        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 3):
            ret_coords.append(-coords[i - 1] + 2 * coords[i] - coords[i + 1])

        ret_coords.append(-coords[-4] + 2 * coords[-3] - coords[-2] - coords[-1])
        ret_coords.append(-coords[-3] + 2 * coords[-2])
        ret_coords.append(-coords[-3] + 2 * coords[-1])
        return tuple(ret_coords)

    class _TypeDOrbitIterator(object):
        '''
        Optimized iterator object that traverses the Weyl group orbit of a given weight.

        Attributes:
            no public attributes
        '''
        def __init__(self, liealg, wt):
            ep_coords = liealg._convert_funds_to_epsilons(liealg.reflect_to_chamber(wt))

            #Check last element to zero/nonzero cases, and store the sign if nonzero
            if ep_coords[-1] < 0:
                self._contains_zero = False
                self._sign = -1
                ep_coords[-1] *= -1
            elif ep_coords[-1] == 0:
                self._contains_zero = True
            else:
                self._contains_zero = False
                self._sign = 1

            #Construct list of items and multiplicities
            self._item_list = [ep_coords[0]]
            cur_item = ep_coords[0]
            rem_list = [0]
            for item in ep_coords:
                if item < cur_item:
                    self._item_list.append(item)
                    rem_list.append(1)
                    cur_item = item
                else:
                    rem_list[-1] += 1

            #Contruct matrix of remaining items, and initial index list
            index_list = []
            rem_mat = [list(rem_list)]
            for i in range(len(ep_coords)):
                j = 0
                while rem_mat[i][j] == 0:
                    j += 1
                index_list.append(j)
                rem_mat.append(list(rem_mat[i]))
                rem_mat[i+1][j] -= 1

            self._index_list = index_list
            self._rem_mat = rem_mat
            self.perms_done = False
            self.liealg = liealg

            # Get first permutation, and construct index iterator, which iterates over subsets of
            # the indices of the non-zero elements of the permutation (except the last element in
            # the case that all elements are nonzero)
            self.cur_perm = self._next_perm()
            if not self._contains_zero:
                non_zero_inds = [i for i in range(len(self.cur_perm) - 1) if self.cur_perm[i] != 0]
            else:
                non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
            self._index_iter = itertools.chain.from_iterable(
                itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds) + 1))

        def __iter__(self):
            return self

        def next(self):
            try:
                neg_inds = self._index_iter.next()
                ep_coords = list(self.cur_perm)
            except StopIteration:
                if self.perms_done:
                    raise StopIteration()
                else:
                    self.cur_perm = self._next_perm()
                    ep_coords = list(self.cur_perm)
                    if not self._contains_zero:
                        non_zero_inds = [i for i in range(len(self.cur_perm) - 1) if self.cur_perm[i] != 0]
                    else:
                        non_zero_inds = [i for i in range(len(self.cur_perm)) if self.cur_perm[i] != 0]
                    self._index_iter = itertools.chain.from_iterable(
                        itertools.combinations(non_zero_inds, r) for r in range(len(non_zero_inds) + 1))
                    neg_inds = self._index_iter.next()

            # Set signs of permutation; last element is negated if there are no zero elements and
            # the number of other negatives is odd
            for i in neg_inds:
                ep_coords[i] = -ep_coords[i]
            if not self._contains_zero:
                ep_coords[-1] *= self._sign * (-1)**len(neg_inds)

            return self.liealg._convert_epsilons_to_funds(ep_coords)

        def _next_perm(self):
            r = len(self._index_list)
            num_items = len(self._item_list)

            # Construct new weight
            ep_coords = []
            for index in self._index_list:
                ep_coords.append(self._item_list[index])

            # Find index to increment
            i = r - 2
            j = 0
            while i >= 0:
                j = self._index_list[i] + 1
                while j < num_items:
                    if self._rem_mat[i][j] > 0: break
                    j += 1
                if j < num_items: break
                i -= 1

            # If we're finished, return the last weight
            if i < 0:
                self.perms_done = True
                return ep_coords

            # Increment indices
            self._index_list[i] = j
            self._rem_mat[i + 1] = list(self._rem_mat[i])
            self._rem_mat[i + 1][j] -= 1
            i += 1
            while i < r:
                j = 0
                while self._rem_mat[i][j] == 0:
                    j += 1
                self._index_list[i] = j
                self._rem_mat[i + 1] = list(self._rem_mat[i])
                self._rem_mat[i + 1][j] -= 1
                i += 1

            return ep_coords

class _Root(tuple):
    """
    This class represents an element of the root lattice of the given simple Lie algebra.
    """

    def __new__(cls, liealg, coords):
        """

        :param liealg:
        :param coords:
        :return:
        """
        return super(_Root, cls).__new__(cls, liealg._convert_roots_to_funds(coords))

    def __init__(self, liealg, coords):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param coords: A list of integers: coordinates for the weight in terms of the basis of simple roots.
        """
        self.liealg = liealg
        self.root_coords = coords

    def get_root_level(self):
        ret_val = 0
        for coord in self.root_coords:
            ret_val = ret_val + coord

        return ret_val
