package lie

import (
	"math/big"
	"sync"

	"github.com/mjschust/cblocks/util"
)

// MutableWeightPoly is a WeightPoly with additional methods to modify the coefficients.
type MutableWeightPoly interface {
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

// NewWeightPolyBuilder constructs a new MutableWeightPoly.
func NewWeightPolyBuilder(rank int) MutableWeightPoly {
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
type PolyProduct func(Weight, Weight) MutableWeightPoly

// MemoReduce reduces the polynomials using memoization.
func (prod PolyProduct) MemoReduce(polys ...WeightPoly) WeightPoly {
	if len(polys) == 0 {
		return nil
	}
	if len(polys) == 1 {
		return polys[0]
	}

	fusionDict := util.NewVectorMap()
	var mutex = &sync.Mutex{}
	var memoProd PolyProduct = func(wt1, wt2 Weight) MutableWeightPoly {
		mutex.Lock()
		submap, present := fusionDict.Get(wt1)
		if present {
			curryMap := submap.(util.VectorMap)
			val, present := curryMap.Get(wt2)
			if present {
				mutex.Unlock()
				return val.(MutableWeightPoly)
			}

			mutex.Unlock()
			rslt := prod(wt1, wt2)
			mutex.Lock()
			curryMap.Put(wt2, rslt)
			mutex.Unlock()
			return rslt
		}

		mutex.Unlock()
		rslt := prod(wt1, wt2)
		mutex.Lock()
		curryMap := util.NewVectorMap()
		curryMap.Put(wt2, rslt)
		fusionDict.Put(wt1, curryMap)
		mutex.Unlock()
		return rslt
	}

	var product = polys[0]
	for i := 1; i < len(polys); i++ {
		product = memoProd.Apply(product, polys[i])
	}

	return product
}

// Reduce applies the product to a list of polynomials
func (prod PolyProduct) Reduce(polys ...WeightPoly) WeightPoly {
	if len(polys) == 0 {
		return nil
	}
	if len(polys) == 1 {
		return polys[0]
	}

	var product = polys[0]
	for i := 1; i < len(polys); i++ {
		product = prod.Apply(product, polys[i])
	}

	return product
}

// Apply the PolyProduct to the given WeightPolys.
func (prod PolyProduct) Apply(poly1, poly2 WeightPoly) WeightPoly {
	type monoExpan struct {
		coeff *big.Int
		poly  WeightPoly
	}
	retPoly := NewWeightPolyBuilder(poly1.Rank())
	poly1Wts := poly1.Weights()
	poly2Wts := poly2.Weights()
	c := make(chan monoExpan)
	for _, wt1 := range poly1Wts {
		mult1 := poly1.Multiplicity(wt1)
		for _, wt2 := range poly2Wts {
			mult2 := poly2.Multiplicity(wt2)
			coeff := big.NewInt(0).Mul(mult1, mult2)
			go func(coeff *big.Int, wt1, wt2 Weight) {
				c <- monoExpan{coeff, prod(wt1, wt2)}
			}(coeff, wt1, wt2)
		}
	}

	coeff := big.NewInt(0)
	for range poly1Wts {
		for range poly2Wts {
			prodRslt := <-c
			summand := prodRslt.poly

			for _, wt := range summand.Weights() {
				mult := summand.Multiplicity(wt)
				coeff.Mul(prodRslt.coeff, mult)
				retPoly.AddMonomial(wt, coeff)
			}
		}
	}
	return retPoly
}
