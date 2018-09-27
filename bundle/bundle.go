package bundle

import (
	"math/big"

	"github.com/mjschust/cblocks/lie"
)

type CBBundle interface {
	Rank() *big.Int
}

func NewSymmetricCBBundle(alg lie.Algebra, wt lie.Weight, ell int, n int) CBBundle {
	wts := make([]lie.Weight, n)
	for i := 0; i < n; i++ {
		wts[i] = wt
	}
	return cbbundleImpl{alg, wts, ell}
}

type cbbundleImpl struct {
	alg lie.Algebra
	wts []lie.Weight
	ell int
}

func (bun cbbundleImpl) Rank() *big.Int {
	alg := bun.alg
	wts := bun.wts
	ell := bun.ell
	fusProd := alg.Fusion(ell, wts[0], wts[1])
	multiFus := alg.Fusion(ell, wts[2:len(wts)]...)

	retVal := big.NewInt(0)
	rslt := big.NewInt(0)
	for _, muStar := range fusProd.Weights() {
		mult1 := fusProd.Multiplicity(muStar)
		mu := alg.Dual(muStar)
		mult2 := multiFus.Multiplicity(mu)
		retVal.Add(retVal, rslt.Mul(mult1, mult2))
	}
	return retVal
}
