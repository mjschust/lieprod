package lie

import "math/big"

// WeightPoly represents a polynomial of weights with integer coefficients.
type WeightPoly interface {
	Rank() int
	Weights() []Weight
	Multiplicity(Weight) *big.Int
}

// Weight represents an integral weight in the weight lattice.
type Weight []int

// Equals compares this weight to another.
func (wt Weight) Equals(other Weight) bool {
	if (wt == nil) != (other == nil) {
		return false
	}

	if len(wt) != len(other) {
		return false
	}

	for i := range wt {
		if wt[i] != other[i] {
			return false
		}
	}

	return true
}

// AddWeights adds the given weights and stores the result in the reciever.
func (wt Weight) AddWeights(wt1, wt2 Weight) {
	for i := range wt1 {
		wt[i] = wt1[i] + wt2[i]
	}
}

// SubWeights subtracts the given weights and stores the result in the reciever.
func (wt Weight) SubWeights(wt1, wt2 Weight) {
	for i := range wt1 {
		wt[i] = wt1[i] - wt2[i]
	}
}

// Rank returns the dimension of the underlying vector.
func (wt Weight) Rank() int {
	return len(wt)
}

// Weights returns a singleton list of this weight.
func (wt Weight) Weights() []Weight {
	return []Weight{wt}
}

// Multiplicity returns one for this weight, zero otherwise.
func (wt Weight) Multiplicity(other Weight) *big.Int {
	if wt.Equals(other) {
		return big.NewInt(1)
	}
	return big.NewInt(0)
}

// Root represents a root in the root lattice.
type Root []int

// Epsilon coordinates are used to perform many root system calculations.
// The basis is type-dependent.
type epCoord []int

func (rslt epCoord) addEpc(wt1, wt2 []int) {
	for i := range wt1 {
		rslt[i] = wt1[i] + wt2[i]
	}
}

func (rslt epCoord) subEpc(wt1, wt2 []int) {
	for i := range wt1 {
		rslt[i] = wt1[i] - wt2[i]
	}
}
