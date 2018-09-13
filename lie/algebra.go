package lie

// Algebra contains type-specific Lie algebra operations.
type Algebra interface {
	DualCoxeter() int
	PositiveRoots() []Root
	Weights(int) []Weight
	Rho() Weight
	KillingForm(Weight, Weight) float64
	IntKillingForm(Weight, Weight) int
	KillingFactor() int
	Level(Weight) int
	Dual(Weight) Weight
	ReflectToChamber(Weight) (Weight, int)
	ConvertRoot(Root) Weight
	NewOrbitIterator(Weight) OrbitIterator
}

// Weight represents a weight.
type Weight []int

// Root represents a root in the root lattice.
type Root []int

// An OrbitIterator supplies an interface to traverse the orbit of a weight.
type OrbitIterator interface {
	Next() Weight
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
	rho := make([]int, 0, alg.rank)
	for i := 0; i < alg.rank; i++ {
		rho = append(rho, 1)
	}

	return rho
}

// KillingForm computes the Killing product of the given weights.
func (alg TypeA) KillingForm(wt1, wt2 Weight) float64 {
	return float64(alg.IntKillingForm(wt1, wt2)) / float64(alg.KillingFactor())
}

// IntKillingForm calculates the Killing product normalized so that the product of integral weights is an integer.
func (alg TypeA) IntKillingForm(wt1, wt2 Weight) int {
	var product int
	epc1 := alg.convertWeight(wt1)
	epc2 := alg.convertWeight(wt2)
	var sum1 int
	var sum2 int

	for i := range epc1 {
		product += epc1[i] * epc2[i]
		sum1 += epc1[i]
		sum2 += epc2[i]
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
func (alg TypeA) Dual(wt Weight) (revWt Weight) {
	revWt = make([]int, len(wt))
	for i := range wt {
		revWt[len(wt)-i-1] = wt[i]
	}
	return
}

// ReflectToChamber reflects the given weight into the dominant chamber and returns
// the result with reflection parity.
func (alg TypeA) ReflectToChamber(wt Weight) (Weight, int) {
	epc := alg.convertWeight(wt)
	epc, parity := insertSort(epc)

	lastCoord := epc[len(epc)-1]
	for i := range epc {
		epc[i] = epc[i] - lastCoord
	}

	return alg.convertEpCoord(epc), parity
}

func insertSort(epc epCoord) (dom epCoord, parity int) {
	dom = make([]int, len(epc))
	copy(dom, epc)
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
	domWt, _ := alg.ReflectToChamber(wt)
	epc := alg.convertWeight(domWt)

	// Construct list of unique coords and multiplicities
	uniqCoords := []int{epc[0]}
	cur := epc[0]
	multiplicities := []int{0}
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
	indices := []int{}
	multMatrix := [][]int{multiplicities}
	for i := range epc {
		j := 0
		for ; multMatrix[i][j] == 0; j++ {
		}
		indices = append(indices, j)
		multMatrix = append(multMatrix, newList(multMatrix[i]))
		multMatrix[i+1][j]--
	}

	return &typeAOrbitIterator{alg, uniqCoords, indices, multMatrix, false}
}

type typeAOrbitIterator struct {
	alg        TypeA
	uniqCoords []int
	indices    []int
	multMatrix [][]int
	done       bool
}

func (iter *typeAOrbitIterator) Next() Weight {
	// Construct new weight
	epc := make([]int, 0, len(iter.indices))
	for _, index := range iter.indices {
		epc = append(epc, iter.uniqCoords[index])
	}
	newWt := iter.alg.convertEpCoord(epc)

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
		return newWt
	}

	// Else, increment indices
	iter.indices[i] = j
	iter.multMatrix[i+1] = newList(iter.multMatrix[i])
	iter.multMatrix[i+1][j]--
	i++
	for ; i < len(iter.indices); i++ {
		j := 0
		for ; iter.multMatrix[i][j] == 0; j++ {
		}
		iter.indices[i] = j
		iter.multMatrix[i+1] = newList(iter.multMatrix[i])
		iter.multMatrix[i+1][j]--
	}

	return newWt
}

func (iter *typeAOrbitIterator) HasNext() bool {
	return !iter.done
}

func newList(list []int) (newCopy []int) {
	newCopy = make([]int, len(list))
	copy(newCopy, list)
	return
}

type epCoord []int

func (alg TypeA) convertWeight(wt Weight) (epc epCoord) {
	epc = make([]int, alg.rank+1)
	var part int
	for i := len(wt) - 1; i >= 0; i-- {
		part += wt[i]
		epc[i] = part
	}
	return
}

func (alg TypeA) convertEpCoord(epc epCoord) (wt Weight) {
	wt = make([]int, alg.rank)
	for i := range wt {
		wt[i] = epc[i] - epc[i+1]
	}
	return
}

// ConvertRoot converts a root into a weight.
func (alg TypeA) ConvertRoot(rt Root) Weight {
	if alg.rank == 1 {
		return []int{2 * rt[0]}
	}

	wt := []int{}
	wt = append(wt, 2*rt[0]-rt[1])
	for i := 1; i < len(rt)-1; i++ {
		wt = append(wt, 2*rt[i]-rt[i+1]-rt[i-1])
	}

	wt = append(wt, 2*rt[len(rt)-1]-rt[len(rt)-2])
	return wt
}
