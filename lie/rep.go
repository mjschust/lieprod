package lie

import (
	"sort"

	"github.com/mjschust/cblocks/util"
)

// ReprDimension returns the dimension of the irreducible representation for the
// given highest weight.
func ReprDimension(alg Algebra, highestWt Weight) int {
	rho := alg.Rho()
	posRoots := alg.PositiveRoots()

	numer := 1
	denom := 1
	for _, root := range posRoots {
		wtForm := alg.ConvertRoot(root)
		a := alg.IntKillingForm(highestWt, wtForm)
		b := alg.IntKillingForm(rho, wtForm)
		numer *= a + b
		denom *= b
	}

	return numer / denom
}

// DominantChar builds the dominant character of the representation of the given highest weight.
func DominantChar(alg Algebra, highestWt Weight) util.VectorMap {
	// Construct root-level map
	posRoots := alg.PositiveRoots()
	rootLevelMap := make(map[int][]Weight)
	for _, root := range posRoots {
		level := 0
		for i := range root {
			level += root[i]
		}

		wtForm := alg.ConvertRoot(root)
		rtList, present := rootLevelMap[level]
		if present {
			rootLevelMap[level] = append(rtList, wtForm)
		} else {
			rootLevelMap[level] = []Weight{wtForm}
		}
	}

	// Construct set of dominant weights
	level := 0
	weightLevelDict := make(map[int]util.VectorMap)
	weightLevelDict[0] = util.NewVectorMap()
	weightLevelDict[0].Put(highestWt, true)
	domWts := util.NewVectorMap()
	domWts.Put(highestWt, true)
	for {
		done := true
		for key := range weightLevelDict {
			if level <= key {
				done = false
				break
			}
		}
		if done {
			break
		}

		_, present := weightLevelDict[level]
		if !present {
			level++
			continue
		}

		for _, wt := range weightLevelDict[level].Keys() {
			for rootLevel, roots := range rootLevelMap {
				for _, root := range roots {
					newWt := subWeights(wt, root)
					if isDominant(newWt) {
						wtSet, present := weightLevelDict[level+rootLevel]
						if present {
							wtSet.Put(newWt, true)
						} else {
							wtSet = util.NewVectorMap()
							wtSet.Put(newWt, true)
							weightLevelDict[level+rootLevel] = wtSet
						}
						domWts.Put(newWt, true)
					}
				}
			}
		}

		level++
	}

	// Build dominant character
	sortedLevels := make([]int, 0, len(weightLevelDict))
	for level := range weightLevelDict {
		sortedLevels = append(sortedLevels, level)
	}
	sort.Slice(sortedLevels, func(i, j int) bool { return sortedLevels[i] < sortedLevels[j] })

	domChar := util.NewVectorMap()
	domChar.Put(highestWt, 1)
	rho := alg.Rho()
	for _, level := range sortedLevels {
		for _, wt := range weightLevelDict[level].Keys() {
			var freudenthalHelper func(wt Weight)
			freudenthalHelper = func(wt Weight) {
				_, present := domChar.Get(wt)
				if present {
					return
				}

				multiplicitySum := 0
				for _, root := range posRoots {
					rootWt := alg.ConvertRoot(root)
					n := 0
					newWeight := wt
					a := alg.IntKillingForm(wt, rootWt)
					b := alg.IntKillingForm(rootWt, rootWt)

					for {
						n++
						newWeight = addWeights(newWeight, rootWt)
						newDomWeight, _ := alg.ReflectToChamber(newWeight)
						_, present := domWts.Get(newDomWeight)
						if !present {
							break
						}

						freudenthalHelper(newDomWeight)
						newWeightMultiplicity, _ := domChar.Get(newDomWeight)
						multiplicitySum += (a + n*b) * newWeightMultiplicity.(int)
					}
				}

				numerator := 2 * multiplicitySum
				tempWt := addWeights(highestWt, rho)
				denominator := alg.IntKillingForm(tempWt, tempWt)
				tempWt = addWeights(wt, rho)
				denominator -= alg.IntKillingForm(tempWt, tempWt)
				domChar.Put(wt, numerator/denominator)
			}
			freudenthalHelper(wt)
		}
	}

	return domChar
}

// Tensor computes the tensor product decomposition of the given representations.
func Tensor(alg Algebra, wt1, wt2 Weight) util.VectorMap {
	if ReprDimension(alg, wt1) < ReprDimension(alg, wt2) {
		wt1, wt2 = wt2, wt1
	}

	rho := alg.Rho()
	domChar := DominantChar(alg, wt2)
	lamRhoSum := addWeights(wt1, rho)
	retDict := util.NewVectorMap()

	for _, domWeight := range domChar.Keys() {
		value, _ := domChar.Get(domWeight)
		domWtMult := value.(int)
		orbitIterator := alg.NewOrbitIterator(domWeight)
		for orbitIterator.HasNext() {
			orbitWeight := orbitIterator.Next()
			newSum := addWeights(lamRhoSum, orbitWeight)
			newDomWeight, parity := alg.ReflectToChamber(newSum)
			newDomWeight = subWeights(newDomWeight, rho)
			if !isDominant(newDomWeight) {
				continue
			}
			value, present := retDict.Get(newDomWeight)
			if present {
				curMult := value.(int)
				retDict.Put(newDomWeight, curMult+parity*domWtMult)
			} else {
				retDict.Put(newDomWeight, parity*domWtMult)
			}
		}
	}

	return retDict
}

func addWeights(wt1, wt2 Weight) Weight {
	var retWt Weight = make([]int, len(wt1))

	for i := range wt1 {
		retWt[i] = wt1[i] + wt2[i]
	}

	return retWt
}

func subWeights(wt1, wt2 Weight) Weight {
	var retWt Weight = make([]int, len(wt1))

	for i := range wt1 {
		retWt[i] = wt1[i] - wt2[i]
	}

	return retWt
}

func isDominant(wt Weight) bool {
	for _, coord := range wt {
		if coord < 0 {
			return false
		}
	}

	return true
}
