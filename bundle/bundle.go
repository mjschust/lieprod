package bundle

import (
	"math/big"

	"github.com/mjschust/cblocks/lie"
)

// CBBundle represents a conformal blocks bundle.
type CBBundle interface {
	Algebra() lie.Algebra
	Weights() []lie.Weight
	Level() int
	Points() int
	Rank() *big.Int
}

// NewCBBundle creates a new CBBundle.
func NewCBBundle(alg lie.Algebra, wts []lie.Weight, ell int) CBBundle {
	n := len(wts)
	newWts := make([]lie.Weight, n)
	for i := 0; i < n; i++ {
		wtCopy := alg.NewWeight()
		copy(wtCopy, wts[i])
		newWts[i] = wtCopy
	}
	return cbbundleImpl{alg, newWts, ell}
}

// NewSymmetricCBBundle creates a CBBundle with the given weight repeated n times.
func NewSymmetricCBBundle(alg lie.Algebra, wt lie.Weight, ell int, n int) CBBundle {
	wts := make([]lie.Weight, n)
	wtCopy := alg.NewWeight()
	copy(wtCopy, wt)
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

func (bun cbbundleImpl) Algebra() lie.Algebra {
	return bun.alg
}

func (bun cbbundleImpl) Weights() []lie.Weight {
	wts := bun.wts
	n := len(wts)
	newWts := make([]lie.Weight, n)
	for i := 0; i < n; i++ {
		wtCopy := bun.alg.NewWeight()
		copy(wtCopy, wts[i])
		newWts[i] = wtCopy
	}
	return newWts
}

func (bun cbbundleImpl) Level() int {
	return bun.ell
}

func (bun cbbundleImpl) Points() int {
	return len(bun.wts)
}

func (bun cbbundleImpl) Rank() *big.Int {
	alg := bun.alg
	wts := bun.wts
	ell := bun.ell
	product := alg.Fusion(ell, wts[1:len(wts)]...)

	return product.Multiplicity(alg.Dual(wts[0]))
}
