package lie

import (
	"github.com/mjschust/cblocks/util"
)

// ReprDimension returns the dimension of the irreducible representation for the
// given highest weight.
func ReprDimension(alg Algebra, highestWt Weight) float64 {
	rho := alg.Rho()
	posRoots := alg.PositiveRoots()

	numer := 1.0
	denom := 1.0
	for _, root := range posRoots {
		wtForm := alg.ConvertRoot(root)
		a := alg.KillingForm(highestWt, wtForm)
		b := alg.KillingForm(rho, wtForm)
		numer *= a + b
		denom *= b
	}

	return numer / denom
}

func DominantWeights(alg Algebra, highestWt Weight) util.VectorMap {
	// Construct root-level map
	posRoots := alg.PositiveRoots()
	rootLevelMap := make(map[float64][]Weight)
	for _, root := range posRoots {
		level := 0.0
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
	level := 0.0
	weightLevelDict := make(map[float64][]Weight)
	weightLevelDict[0] = []Weight{highestWt}
	domWts := util.NewVectorMap()
	domWts.Put(highestWt, true)
	for true {
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

		for wt := range weightLevelDict[level] {
			for rootLevel, roots := range rootLevelMap {
				for root := range roots {
					// TODO
				}
			}
		}
	}

	return domWts
}
