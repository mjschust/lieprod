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

// A WeightProduct defines a product on weights with polynomial output.
type WeightProduct func(Weight, Weight) MutableWeightPoly

// A PolyProduct defines a product on weight polynomials.
type PolyProduct interface {
	Apply(WeightPoly, WeightPoly) WeightPoly
	Reduce(...WeightPoly) WeightPoly
}

// NewProduct constructs a poly product without memoization.
func NewProduct(prod WeightProduct) PolyProduct {
	return polyProductImpl{plainKernel{prod}}
}

// NewMemoizedProduct constructs a poly product with memoization.
func NewMemoizedProduct(prod WeightProduct) PolyProduct {
	return polyProductImpl{&memoizedKernel{prod, util.NewVectorMap(), sync.Mutex{}}}
}

type polyProductImpl struct {
	productKernel
}

func (app polyProductImpl) Apply(poly1, poly2 WeightPoly) WeightPoly {
	type monoExpan struct {
		coeff   *big.Int
		promise polyPromise
	}
	retPoly := NewWeightPolyBuilder(poly1.Rank())
	poly1Wts := poly1.Weights()
	poly2Wts := poly2.Weights()
	expanList := make([]monoExpan, 0, len(poly1Wts)*len(poly2Wts))
	for _, wt1 := range poly1Wts {
		mult1 := poly1.Multiplicity(wt1)
		for _, wt2 := range poly2Wts {
			mult2 := poly2.Multiplicity(wt2)
			coeff := big.NewInt(0).Mul(mult1, mult2)
			promise := app.asynchApply(wt1, wt2)
			expanList = append(expanList, monoExpan{coeff, promise})
		}
	}

	rslt := big.NewInt(0)
	for _, exp := range expanList {
		coeff := exp.coeff
		summand := exp.promise()

		for _, wt := range summand.Weights() {
			mult := summand.Multiplicity(wt)
			rslt.Mul(coeff, mult)
			retPoly.AddMonomial(wt, rslt)
		}
	}

	return retPoly
}

func (app polyProductImpl) Reduce(polys ...WeightPoly) WeightPoly {
	if len(polys) == 0 {
		return nil
	}
	if len(polys) == 1 {
		return polys[0]
	}

	var product = polys[0]
	for i := 1; i < len(polys); i++ {
		product = app.Apply(product, polys[i])
	}

	return product
}

type productKernel interface {
	asynchApply(Weight, Weight) polyPromise
}

type polyPromise func() WeightPoly

type plainKernel struct {
	prod WeightProduct
}

func (knl plainKernel) asynchApply(wt1, wt2 Weight) polyPromise {
	c := make(chan WeightPoly)
	go func(wt1, wt2 Weight, c chan WeightPoly) {
		c <- knl.prod(wt1, wt2)
	}(wt1, wt2, c)

	return func() WeightPoly {
		rslt := <-c
		return rslt
	}
}

type memoizedKernel struct {
	prod     WeightProduct
	rsltDict util.VectorMap
	sync.Mutex
}

func (knl *memoizedKernel) asynchApply(wt1, wt2 Weight) polyPromise {
	knl.Lock()
	submap, present := knl.rsltDict.Get(wt1)
	if present {
		curryMap := submap.(util.VectorMap)
		val, present := curryMap.Get(wt2)
		if present {
			knl.Unlock()
			return func() WeightPoly {
				return val.(WeightPoly)
			}
		}

		knl.Unlock()
		c := make(chan WeightPoly)
		go func(wt1, wt2 Weight, c chan WeightPoly) {
			c <- knl.prod(wt1, wt2)
		}(wt1, wt2, c)

		return func() WeightPoly {
			rslt := <-c
			knl.Lock()
			curryMap.Put(wt2, rslt)
			knl.Unlock()
			return rslt
		}
	}

	knl.Unlock()
	c := make(chan WeightPoly)
	go func(wt1, wt2 Weight, c chan WeightPoly) {
		c <- knl.prod(wt1, wt2)
	}(wt1, wt2, c)

	return func() WeightPoly {
		rslt := <-c
		knl.Lock()
		submap, present := knl.rsltDict.Get(wt1)
		if present {
			curryMap := submap.(util.VectorMap)
			curryMap.Put(wt2, rslt)
			knl.Unlock()
			return rslt
		}

		curryMap := util.NewVectorMap()
		curryMap.Put(wt2, rslt)
		knl.rsltDict.Put(wt1, curryMap)
		knl.Unlock()
		return rslt
	}
}
