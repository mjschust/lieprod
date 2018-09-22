package lie

import (
	"math/big"
	"sort"

	"github.com/mjschust/cblocks/util"
)

// An Algebra supplies representation-theoretic methods acting on Weights.
type Algebra interface {
	ReprDimension(Weight) *big.Int
	DominantChar(Weight) WeightPoly
	Tensor(Weight, Weight) WeightPoly
}

type algebraImpl struct {
	RootSystem
}

// NewAlgebra constructs and returns the lie algebra associated to the given root system.
func NewAlgebra(rtsys RootSystem) Algebra {
	return algebraImpl{rtsys}
}

// ReprDimension returns the dimension of the irreducible representation for the
// given highest weight.
func (alg algebraImpl) ReprDimension(highestWt Weight) *big.Int {
	rho := alg.Rho()
	posRoots := alg.PositiveRoots()

	numer := big.NewInt(1)
	denom := big.NewInt(1)
	a := big.NewInt(0)
	b := big.NewInt(0)
	rslt := big.NewInt(0)
	wtForm := alg.NewWeight()
	for _, root := range posRoots {
		alg.convertRoot(root, wtForm)
		a.SetInt64(int64(alg.IntKillingForm(highestWt, wtForm)))
		b.SetInt64(int64(alg.IntKillingForm(rho, wtForm)))
		numer.Mul(numer, rslt.Add(a, b))
		denom.Mul(denom, b)
	}

	return rslt.Div(numer, denom)
}

// DominantChar builds the dominant character of the representation of the given highest weight.
func (alg algebraImpl) DominantChar(highestWt Weight) WeightPoly {
	// Construct root-level map
	posRoots := alg.PositiveRoots()
	rootLevelMap := make(map[int][]Weight)
	for _, root := range posRoots {
		level := 0
		for i := range root {
			level += root[i]
		}

		wtForm := alg.NewWeight()
		alg.convertRoot(root, wtForm)
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
	domChar := NewWeightPolyBuilder(alg.Rank())
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

		newWt := alg.NewWeight()
		for _, wt := range weightLevelDict[level].Keys() {
			for rootLevel, roots := range rootLevelMap {
				for _, root := range roots {
					newWt.SubWeights(wt, root)
					if isDominant(newWt) {
						polyWeight := domChar.addWeight(newWt)
						wtSet, present := weightLevelDict[level+rootLevel]
						if present {
							wtSet.Put(polyWeight, true)
						} else {
							wtSet = util.NewVectorMap()
							wtSet.Put(polyWeight, true)
							weightLevelDict[level+rootLevel] = wtSet
						}
						domChar.SetMultiplicity(polyWeight, big.NewInt(-1))
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

	one := big.NewInt(1)
	domChar.SetMultiplicity(highestWt, one)
	rho := alg.Rho()
	for _, level := range sortedLevels {
		for _, wt := range weightLevelDict[level].Keys() {
			var freudenthalHelper func(wt Weight)
			freudenthalHelper = func(wt Weight) {
				mult := domChar.Multiplicity(wt)
				if mult.Sign() > 0 {
					return
				}

				multiplicitySum := big.NewInt(0)
				a := big.NewInt(0)
				b := big.NewInt(0)
				n := big.NewInt(0)
				rslt := big.NewInt(0)
				shiftedWeight := alg.NewWeight()
				newDomWeight := alg.NewWeight()
				workingEpc := alg.newEpc()
				for _, roots := range rootLevelMap {
					for _, rootWt := range roots {
						n.SetInt64(0)
						copy(shiftedWeight, wt)
						a.SetInt64(int64(alg.IntKillingForm(wt, rootWt)))
						b.SetInt64(int64(alg.IntKillingForm(rootWt, rootWt)))

						for {
							n.Add(n, one)
							shiftedWeight.AddWeights(shiftedWeight, rootWt)
							alg.convertWeightToEpc(shiftedWeight, workingEpc)
							alg.reflectEpcToChamber(workingEpc)
							alg.convertEpCoord(workingEpc, newDomWeight)
							if domChar.Multiplicity(newDomWeight).Sign() == 0 {
								break
							}

							freudenthalHelper(newDomWeight)
							newWeightMultiplicity := domChar.Multiplicity(newDomWeight)
							rslt.Add(a, rslt.Mul(n, b))
							rslt.Mul(rslt, newWeightMultiplicity)
							multiplicitySum.Add(multiplicitySum, rslt)
						}
					}
				}

				numerator := big.NewInt(2)
				numerator.Mul(numerator, multiplicitySum)
				shiftedWeight.AddWeights(highestWt, rho)
				denominator := big.NewInt(int64(
					alg.IntKillingForm(shiftedWeight, shiftedWeight)))
				shiftedWeight.AddWeights(wt, rho)
				rslt.SetInt64(int64(alg.IntKillingForm(shiftedWeight, shiftedWeight)))
				denominator.Sub(denominator, rslt)
				domChar.SetMultiplicity(wt, rslt.Div(numerator, denominator))
			}
			freudenthalHelper(wt)
		}
	}

	return domChar
}

// Tensor computes the tensor product decomposition of the given representations.
func (alg algebraImpl) Tensor(wt1, wt2 Weight) WeightPoly {
	if alg.ReprDimension(wt1).Cmp(alg.ReprDimension(wt2)) < 0 {
		wt1, wt2 = wt2, wt1
	}

	// Construct constant weights
	rho := alg.newEpc()
	alg.convertWeightToEpc(alg.Rho(), rho)
	domChar := alg.DominantChar(wt2)
	lamRhoSumWt := alg.NewWeight()
	lamRhoSumWt.AddWeights(wt1, alg.Rho())
	lamRhoSum := alg.newEpc()
	alg.convertWeightToEpc(lamRhoSumWt, lamRhoSum)

	// Construct return map
	retPoly := NewWeightPolyBuilder(alg.Rank())
	var epc = alg.newEpc()
	var orbitEpc = alg.newEpc()
	domWeight := alg.NewWeight()
	rslt := big.NewInt(0)
	for _, charWeight := range domChar.Weights() {
		domWtMult := domChar.Multiplicity(charWeight)
		alg.convertWeightToEpc(charWeight, orbitEpc)
		done := false
		for ; !done; done = alg.nextOrbitEpc(orbitEpc) {
			// Shifted reflection into dominant chamber
			epc.addEpc(lamRhoSum, orbitEpc)
			parity := alg.reflectEpcToChamber(epc)
			epc.subEpc(epc, rho)

			// Check if dominant
			alg.convertEpCoord(epc, domWeight)
			if !isDominant(domWeight) {
				continue
			}

			// Set new multiplicity
			rslt.SetInt64(int64(parity))
			rslt.Mul(rslt, domWtMult)
			retPoly.AddMultiplicity(domWeight, rslt)
		}
	}

	return retPoly
}

func isDominant(wt Weight) bool {
	for _, coord := range wt {
		if coord < 0 {
			return false
		}
	}

	return true
}
