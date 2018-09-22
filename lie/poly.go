package lie

import (
	"math/big"

	"github.com/mjschust/cblocks/util"
)

// WeightPolyBuilder is a WeightPoly with additional methods to modify the coefficients.
type WeightPolyBuilder interface {
	WeightPoly
	addWeight(wt Weight) Weight
	SetMultiplicity(Weight, *big.Int)
	AddMultiplicity(Weight, *big.Int)
}

type hashPolyBuilder struct {
	rank int
	vmap util.VectorMap
}

// NewWeightPolyBuilder constructs a new WeightPolyBuilder.
func NewWeightPolyBuilder(rank int) WeightPolyBuilder {
	return hashPolyBuilder{rank, util.NewVectorMap()}
}

func (poly hashPolyBuilder) Weights() []Weight {
	keys := poly.vmap.Keys()
	retSlc := make([]Weight, len(keys))
	for i, key := range keys {
		retSlc[i] = key
	}
	return retSlc
}

func (poly hashPolyBuilder) Multiplicity(wt Weight) *big.Int {
	val, present := poly.vmap.Get(wt)
	if present {
		return val.(*big.Int)
	}
	return big.NewInt(0)
}

func (poly hashPolyBuilder) addWeight(wt Weight) Weight {
	_, present := poly.vmap.Get(wt)
	if !present {
		newWt := make([]int, poly.rank)
		copy(newWt, wt)
		poly.vmap.Put(newWt, big.NewInt(0))
		return newWt
	}
	return wt
}

func (poly hashPolyBuilder) SetMultiplicity(wt Weight, val *big.Int) {
	curVal, present := poly.vmap.Get(wt)
	if !present {
		poly.addWeight(wt)
		poly.SetMultiplicity(wt, val)
		return
	}

	poly.vmap.Put(wt, curVal.(*big.Int).Set(val))
}

func (poly hashPolyBuilder) AddMultiplicity(wt Weight, val *big.Int) {
	curVal, present := poly.vmap.Get(wt)
	if !present {
		poly.addWeight(wt)
		poly.AddMultiplicity(wt, val)
		return
	}

	poly.vmap.Put(wt, curVal.(*big.Int).Add(curVal.(*big.Int), val))
}
