package lie

// Algebra contains type-specific Lie algebra operations.
type Algebra interface {
	DualCoxeter() int
	NewWeight() Weight
	PositiveRoots() []Root
	Weights(int) []Weight
	Rho() Weight
	KillingForm(Weight, Weight) float64
	IntKillingForm(Weight, Weight) int
	KillingFactor() int
	Level(Weight) int
	dual(Weight, Weight)
	reflectToChamber(Weight, Weight) int
	convertRoot(Root, Weight)
	NewOrbitIterator(Weight) OrbitIterator
}

// Weight represents a weight.
type Weight []int

// Dual computes the dual of the given weight, with the result stored in the reciever.
// WARNING: not safe for in-place usage.
func (rslt Weight) Dual(alg Algebra, wt Weight) {
	alg.dual(wt, rslt)
}

// ConvertRoot converts the given root into a weight, with the result stored in the reciever.
func (rslt Weight) ConvertRoot(alg Algebra, rt Root) {
	alg.convertRoot(rt, rslt)
}

// ReflectToChamber reflects wt into the dominant chamber, storing the result in the reciever
// and returning the parity of the reflection.
func (rslt Weight) ReflectToChamber(alg Algebra, wt Weight) int {
	return alg.reflectToChamber(wt, rslt)
}

// AddWeights adds the given weights and stores the result in the reciever.
func (rslt Weight) AddWeights(wt1, wt2 Weight) {
	for i := range wt1 {
		rslt[i] = wt1[i] + wt2[i]
	}
}

// SubWeights subtracts the given weights and stores the result in the reciever.
func (rslt Weight) SubWeights(wt1, wt2 Weight) {
	for i := range wt1 {
		rslt[i] = wt1[i] - wt2[i]
	}
}

// Root represents a root in the root lattice.
type Root []int

// An OrbitIterator supplies an interface to traverse the orbit of a weight.
type OrbitIterator interface {
	Next(Weight) Weight
	HasNext() bool
}

// TypeA represents the Lie algebra of type A with the specified rank.
type TypeA struct {
	rank int
}

// DualCoxeter computes the dual Coxeter number of the Lie algebra.
func (alg TypeA) DualCoxeter() int {
	return alg.rank + 1
}

// NewWeight creates a new zero weight.
func (alg TypeA) NewWeight() Weight {
	return make([]int, alg.rank, alg.rank+1)
}

// PositiveRoots builds a list of all positive roots of the Lie algebra.
func (alg TypeA) PositiveRoots() []Root {
	retList := make([]Root, 0, alg.rank*(alg.rank+1)/2)
	root := make([]int, alg.rank)

	for i := range root {
		for j := i; j < len(root); j++ {
			root[j] = 1

			var next Root = make([]int, alg.rank)
			copy(next, root)
			retList = append(retList, next)
		}

		for j := i; j < len(root); j++ {
			root[j] = 0
		}
	}

	return retList
}

// Weights returns a slice of all weights with level at most the given int.
func (alg TypeA) Weights(level int) []Weight {
	var weightsHelper func(rank int) []Weight
	weightsHelper = func(rank int) []Weight {
		retList := make([]Weight, 0)
		if rank == 1 {
			for i := 0; i <= level; i++ {
				retList = append(retList, Weight{i})
			}
		} else {
			rMinusOneList := weightsHelper(rank - 1)
			for _, wt := range rMinusOneList {
				wtLevel := 0
				for _, coord := range wt {
					wtLevel += coord
				}
				for i := 0; i <= level-wtLevel; i++ {
					newWt := make([]int, len(wt)+1)
					copy(newWt, wt)
					newWt[len(wt)] = i
					retList = append(retList, newWt)
				}
			}
		}

		return retList
	}

	return weightsHelper(alg.rank)
}

// Rho returns one-half the sum of the positive roots of the algebra.
func (alg TypeA) Rho() Weight {
	rho := alg.NewWeight()
	for i := 0; i < alg.rank; i++ {
		rho[i] = 1
	}

	return rho
}

// KillingForm computes the Killing product of the given weights.
func (alg TypeA) KillingForm(wt1, wt2 Weight) float64 {
	return float64(alg.IntKillingForm(wt1, wt2)) / float64(alg.KillingFactor())
}

// IntKillingForm calculates the Killing product normalized so that the product of integral weights is an integer.
func (alg TypeA) IntKillingForm(wt1, wt2 Weight) int {
	var part1, part2, product, sum1, sum2 int

	for i := len(wt1) - 1; i >= 0; i-- {
		part1 += wt1[i]
		part2 += wt2[i]
		product += part1 * part2
		sum1 += part1
		sum2 += part2
	}

	return (alg.rank+1)*product - sum1*sum2
}

// KillingFactor returns IntKillingForm/KillingForm.
func (alg TypeA) KillingFactor() int {
	return alg.rank + 1
}

// Level computes the 'level' of the given weight, i.e. its product with the highest root.
func (alg TypeA) Level(wt Weight) (lv int) {
	for i := range wt {
		lv += wt[i]
	}
	return
}

// Dual computes the highest weight of the dual repr. of corresponding to the given weight.
func (alg TypeA) dual(wt Weight, rslt Weight) {
	for i := range wt {
		rslt[len(wt)-i-1] = wt[i]
	}
}

// ReflectToChamber reflects the given weight into the dominant chamber and returns
// the result with reflection parity.
func (alg TypeA) reflectToChamber(wt Weight, rslt Weight) int {
	var epc []int
	if cap(rslt) > len(rslt) {
		epc = append(rslt, 0)
		epc = alg.convertWeightToEpc(wt, epc)
	} else {
		epc = make([]int, alg.rank+1)
		epc = alg.convertWeightToEpc(wt, epc)
	}
	epc, parity := insertSort(epc)

	lastCoord := epc[len(epc)-1]
	for i := range epc {
		epc[i] = epc[i] - lastCoord
	}

	alg.convertEpCoord(epc, rslt)
	return parity
}

func insertSort(epc []int) (dom []int, parity int) {
	dom = epc[:]
	parity = 1

	for i := range dom {
		for j := i; j > 0 && dom[j-1] < dom[j]; j-- {
			dom[j-1], dom[j] = dom[j], dom[j-1]
			parity *= -1
		}
	}
	return
}

// NewOrbitIterator creates a new OrbitIterator for the given weight.
func (alg TypeA) NewOrbitIterator(wt Weight) OrbitIterator {
	domWt := alg.NewWeight()
	domWt.ReflectToChamber(alg, wt)
	epc := alg.convertWeightToEpc(domWt, append(domWt, 0))

	// Construct list of unique coords and multiplicities
	uniqCoords := make([]int, 0, len(epc))
	uniqCoords = append(uniqCoords, epc[0])
	cur := epc[0]
	multiplicities := make([]int, 0, len(epc))
	multiplicities = append(multiplicities, 0)
	for _, coord := range epc {
		if coord < cur {
			uniqCoords = append(uniqCoords, coord)
			multiplicities = append(multiplicities, 1)
			cur = coord
		} else {
			multiplicities[len(multiplicities)-1]++
		}
	}

	// Initialize multiplicity matrix and indices
	indices := make([]int, 0, len(epc))
	multMatrix := make([][]int, 0, len(epc)+1)
	multMatrix = append(multMatrix, multiplicities)
	for i := range epc {
		j := 0
		for ; multMatrix[i][j] == 0; j++ {
		}
		indices = append(indices, j)

		newCopy := make([]int, len(multMatrix[i]))
		copy(newCopy, multMatrix[i])
		multMatrix = append(multMatrix, newCopy)
		multMatrix[i+1][j]--
	}

	return &typeAOrbitIterator{alg, uniqCoords, indices, multMatrix, epc, false}
}

type typeAOrbitIterator struct {
	alg        TypeA
	uniqCoords []int
	indices    []int
	multMatrix [][]int
	epc        []int
	done       bool
}

func (iter *typeAOrbitIterator) Next(rslt Weight) Weight {
	// Construct new weight
	for i, index := range iter.indices {
		iter.epc[i] = iter.uniqCoords[index]
	}
	iter.alg.convertEpCoord(iter.epc, rslt)

	// Find index to increment
	i := len(iter.indices) - 2
	j := 0
	for ; i >= 0; i-- {
		j = iter.indices[i] + 1
		for ; j < len(iter.uniqCoords); j++ {
			if iter.multMatrix[i][j] > 0 {
				break
			}
		}
		if j < len(iter.uniqCoords) {
			break
		}
	}

	// If we're finished, return the last weight
	if i < 0 {
		iter.done = true
		return rslt
	}

	// Else, increment indices
	iter.indices[i] = j
	copy(iter.multMatrix[i+1], iter.multMatrix[i])
	iter.multMatrix[i+1][j]--
	i++
	for ; i < len(iter.indices); i++ {
		j := 0
		for ; iter.multMatrix[i][j] == 0; j++ {
		}
		iter.indices[i] = j
		copy(iter.multMatrix[i+1], iter.multMatrix[i])
		iter.multMatrix[i+1][j]--
	}

	return rslt
}

func (iter *typeAOrbitIterator) HasNext() bool {
	return !iter.done
}

func (alg TypeA) convertWeightToEpc(wt Weight, epc []int) []int {
	var part int
	for i := len(wt) - 1; i >= 0; i-- {
		part += wt[i]
		epc[i] = part
	}

	return epc
}

func (alg TypeA) convertEpCoord(epc []int, retVal Weight) {
	for i := range retVal {
		retVal[i] = epc[i] - epc[i+1]
	}
}

// ConvertRoot converts a root into a weight.
func (alg TypeA) convertRoot(rt Root, rslt Weight) {
	if alg.rank == 1 {
		rslt[0] = 2 * rt[0]
		return
	}

	rslt[0] = 2*rt[0] - rt[1]
	for i := 1; i < len(rt)-1; i++ {
		rslt[i] = 2*rt[i] - rt[i+1] - rt[i-1]
	}

	rslt[len(rt)-1] = 2*rt[len(rt)-1] - rt[len(rt)-2]
}
