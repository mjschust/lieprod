package lie

import (
	"math/big"

	"github.com/mjschust/cblocks/util"
)

// WeightPolyBuilder is a WeightPoly with additional methods to modify the coefficients.
type WeightPolyBuilder interface {
	WeightPoly
	addWeight(wt Weight) Weight
	SetMonomial(Weight, *big.Int)
	AddMonomial(Weight, *big.Int)
	Add(WeightPoly)
	Mult(*big.Int)
}

type hashPolyBuilder struct {
	rank int
	vmap util.VectorMap
}

// NewWeightPolyBuilder constructs a new WeightPolyBuilder.
func NewWeightPolyBuilder(rank int) WeightPolyBuilder {
	return hashPolyBuilder{rank, util.NewVectorMap()}
}

func (poly hashPolyBuilder) Rank() int {
	return poly.rank
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

func (poly hashPolyBuilder) SetMonomial(wt Weight, val *big.Int) {
	curVal, present := poly.vmap.Get(wt)
	if !present {
		poly.addWeight(wt)
		poly.SetMonomial(wt, val)
		return
	}

	poly.vmap.Put(wt, curVal.(*big.Int).Set(val))
}

func (poly hashPolyBuilder) AddMonomial(wt Weight, val *big.Int) {
	curVal, present := poly.vmap.Get(wt)
	if !present {
		poly.addWeight(wt)
		poly.AddMonomial(wt, val)
		return
	}

	poly.vmap.Put(wt, curVal.(*big.Int).Add(curVal.(*big.Int), val))
}

func (poly hashPolyBuilder) Add(poly2 WeightPoly) {
	for _, wt := range poly2.Weights() {
		mult := poly2.Multiplicity(wt)
		poly.AddMonomial(wt, mult)
	}
}

func (poly hashPolyBuilder) Mult(val *big.Int) {
	for _, wt := range poly.Weights() {
		curVal, _ := poly.vmap.Get(wt)
		poly.vmap.Put(wt, curVal.(*big.Int).Mul(curVal.(*big.Int), val))
	}
}

// A PolyProduct defines a product on WeightPolys.
type PolyProduct func(Weight, Weight) WeightPoly

// Apply the PolyProduct to the given WeightPolys.
func (prod PolyProduct) Apply(poly1, poly2 WeightPoly) WeightPoly {
	retPoly := NewWeightPolyBuilder(poly1.Rank())
	for _, wt1 := range poly1.Weights() {
		mult1 := poly1.Multiplicity(wt1)
		for _, wt2 := range poly2.Weights() {
			mult2 := poly2.Multiplicity(wt2)
			summand := prod(wt1, wt2)
			summand.(WeightPolyBuilder).Mult(mult1)
			summand.(WeightPolyBuilder).Mult(mult2)
			retPoly.Add(summand)
		}
	}
	return retPoly
}
