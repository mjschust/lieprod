package lie

import (
	"math/big"

	"github.com/mjschust/cblocks/util"
)

var zero = big.NewInt(0)

// WeightPolyBuilder is a WeightPoly with additional methods to modify the coefficients.
type WeightPolyBuilder interface {
	WeightPoly
	SetMultiplicity(Weight, *big.Int)
	AddMultiplicity(Weight, *big.Int)
}

type hashPolyBuilder struct {
	rtsys RootSystem
	vmap  util.VectorMap
}

// NewWeightPolyBuilder constructs a new WeightPolyBuilder.
func NewWeightPolyBuilder(rtsys RootSystem) WeightPolyBuilder {
	return hashPolyBuilder{rtsys, util.NewVectorMap()}
}

func (poly hashPolyBuilder) GetWeights() []Weight {
	keys := poly.vmap.Keys()
	retSlc := make([]Weight, len(keys))
	for i, key := range keys {
		retSlc[i] = key
	}
	return retSlc
}

func (poly hashPolyBuilder) GetMultiplicity(wt Weight) *big.Int {
	val, present := poly.vmap.Get(wt)
	if present {
		return val.(*big.Int)
	}
	return zero
}

func (poly hashPolyBuilder) addWeight(wt Weight) {
	_, present := poly.vmap.Get(wt)
	if !present {
		newWt := poly.rtsys.NewWeight()
		copy(newWt, wt)
		poly.vmap.Put(newWt, big.NewInt(0))
	}
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
