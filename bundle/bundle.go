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

// SymCBBundle represents a symmetric conformal blocks bundle.
type SymCBBundle interface {
	CBBundle
	SymmetrizedDivisor() []*big.Rat
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
func NewSymmetricCBBundle(alg lie.Algebra, wt lie.Weight, ell int, n int) SymCBBundle {
	wts := make([]lie.Weight, n)
	wtCopy := alg.NewWeight()
	copy(wtCopy, wt)
	for i := 0; i < n; i++ {
		wts[i] = wt
	}
	return symCbbundleImpl{cbbundleImpl{alg, wts, ell}}
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

type symCbbundleImpl struct {
	cbbundleImpl
}

func (bun symCbbundleImpl) SymmetrizedDivisor() []*big.Rat {
	// Only works for symmetric bundles...need to refactor
	alg := bun.alg
	wts := make([]lie.WeightPoly, len(bun.wts))
	for i := range bun.wts {
		wts[i] = bun.wts[i]
	}
	n := len(bun.wts)
	rho := alg.Rho()
	prod := alg.FusionProduct(bun.ell)
	rk := prod.Reduce(wts[1:n]...).Multiplicity(alg.Dual(bun.wts[0]))

	// Calculate the starting point of each coord
	rslt := big.NewInt(0)
	bun.casimirScalar(wts[0].(lie.Weight), rho, rslt)
	rslt.Mul(rslt, rk)
	denom := big.NewInt(int64(n - 1))
	baseSummand := big.NewRat(0, 1)
	baseSummand.SetFrac(rslt, denom)

	// Prepare the common denominator for the coords
	sumDenom := big.NewRat(0, 1)
	denom.SetInt64(int64(2 * (bun.ell + alg.DualCoxeter()) * alg.KillingFactor()))
	sumDenom.SetInt(denom)

	// Compute the output vector
	a := big.NewRat(0, 1)
	retVec := make([]*big.Rat, n/2-1)
	for i := 2; i < n/2+1; i++ {
		// Compute positive summand
		summand := big.NewRat(0, 1)
		summand.Set(baseSummand)
		a.SetFrac64(int64(i*(n-i)), 1)
		summand.Mul(summand, a)

		// Compute "weighted" factorization calculation
		poly1 := prod.Reduce(wts[0:i]...)
		poly2 := prod.Reduce(wts[i:n]...)
		wfSum := big.NewInt(0)
		for _, mustar := range poly1.Weights() {
			mu := alg.Dual(mustar)
			bun.casimirScalar(mu, rho, rslt)
			rslt.Mul(rslt, poly1.Multiplicity(mustar))
			rslt.Mul(rslt, poly2.Multiplicity(mu))
			wfSum.Add(wfSum, rslt)
		}
		a.SetInt(wfSum)

		summand.Sub(summand, a)
		summand.Quo(summand, sumDenom)
		retVec[i-2] = summand
	}

	return retVec
}

func (bun cbbundleImpl) casimirScalar(wt, rho lie.Weight, rslt *big.Int) {
	rslt.SetInt64(int64(
		bun.alg.IntKillingForm(wt, wt) + 2*bun.alg.IntKillingForm(wt, rho)))
}
