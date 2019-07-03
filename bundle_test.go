package cblocks

import (
	"math/big"
	"testing"

	"github.com/mjschust/cblocks/lie"
)

func TestSymmetricL1Rank(t *testing.T) {
	cases := []struct {
		rtsys lie.RootSystem
		wt    lie.Weight
		ell   int
		n     int
		want  *big.Int
	}{
		// n=3
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 0}, 1, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{1, 0, 0, 0}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 1, 0, 0}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 1, 0}, 1, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 1}, 1, 3, big.NewInt(0)},

		// n=4
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 1, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 1, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 0}, 1, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{1, 0, 0, 0}, 1, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 1, 0, 0}, 1, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 1, 0}, 1, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 1}, 1, 4, big.NewInt(0)},

		// n=5
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 1, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{1, 0, 0, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 1, 0, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 1, 0}, 1, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(4), lie.Weight{0, 0, 0, 1}, 1, 5, big.NewInt(1)},
	}

	for _, c := range cases {
		alg := lie.NewAlgebra(c.rtsys)
		bun := NewSymmetricCBBundle(alg, c.wt, c.ell, c.n)
		got := bun.Rank()
		if got.Cmp(c.want) != 0 {
			t.Errorf("For type A CBBundle(wt: %v, ell: %v, n: %v): Rank() = %v, want %v",
				c.wt, c.ell, c.n, got, c.want)
		}
	}
}

func TestSymmetricL2Rank(t *testing.T) {
	cases := []struct {
		rtsys lie.RootSystem
		wt    lie.Weight
		ell   int
		n     int
		want  *big.Int
	}{
		// n=3
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), lie.Weight{2}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 2}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 0}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 1}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 2, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 0}, 2, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 2, 3, big.NewInt(0)},

		// n=4
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 2, 4, big.NewInt(2)},
		{lie.NewTypeARootSystem(1), lie.Weight{2}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 2, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 2}, 2, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 2, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 2, 4, big.NewInt(2)},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 0}, 2, 4, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 2, 4, big.NewInt(3)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 1}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 2, 4, big.NewInt(3)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 0}, 2, 4, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 2, 4, big.NewInt(1)},

		// n=5
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 2, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), lie.Weight{2}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 2, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 2}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 2, 5, big.NewInt(3)},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 2, 5, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 1}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 2, 5, big.NewInt(5)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 0}, 2, 5, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 2, 5, big.NewInt(0)},
	}

	for _, c := range cases {
		alg := lie.NewAlgebra(c.rtsys)
		bun := NewSymmetricCBBundle(alg, c.wt, c.ell, c.n)
		got := bun.Rank()
		if got.Cmp(c.want) != 0 {
			t.Errorf("For type A CBBundle(wt: %v, ell: %v, n: %v): Rank() = %v, want %v",
				c.wt, c.ell, c.n, got, c.want)
		}
	}
}

func TestSymmetricL3Rank(t *testing.T) {
	cases := []struct {
		rtsys lie.RootSystem
		wt    lie.Weight
		ell   int
		n     int
		want  *big.Int
	}{
		// n=3
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), lie.Weight{2}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), lie.Weight{3}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 2}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 3}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 3, 3, big.NewInt(2)},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 2}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 1}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), lie.Weight{3, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 3}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 2}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 3, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 3, 3, big.NewInt(2)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 2}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 2, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 1}, 3, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 1, 0}, 3, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), lie.Weight{3, 0, 0}, 3, 3, big.NewInt(0)},
	}

	for _, c := range cases {
		alg := lie.NewAlgebra(c.rtsys)
		bun := NewSymmetricCBBundle(alg, c.wt, c.ell, c.n)
		got := bun.Rank()
		if got.Cmp(c.want) != 0 {
			t.Errorf("For type A CBBundle(wt: %v, ell: %v, n: %v): Rank() = %v, want %v",
				c.wt, c.ell, c.n, got, c.want)
		}
	}
}

func TestL3RandomRank(t *testing.T) {
	cases := []struct {
		rtsys lie.RootSystem
		wts   []lie.Weight
		ell   int
		want  *big.Int
	}{
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{0}, lie.Weight{1}, lie.Weight{1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{2}, lie.Weight{1}, lie.Weight{2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{3}, lie.Weight{1}, lie.Weight{3}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{0}, lie.Weight{0}, lie.Weight{0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{3}, lie.Weight{0}, lie.Weight{0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{2}, lie.Weight{0}, lie.Weight{2}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{3}, lie.Weight{0}, lie.Weight{1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{2}, lie.Weight{3}, lie.Weight{0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{3}, lie.Weight{3}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{3}, lie.Weight{3}, lie.Weight{1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 2}, lie.Weight{2, 1}, lie.Weight{2, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 0}, lie.Weight{2, 1}, lie.Weight{0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{0, 3}, lie.Weight{0, 3}, lie.Weight{0, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{3, 0}, lie.Weight{2, 0}, lie.Weight{2, 1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 1}, lie.Weight{1, 1}, lie.Weight{3, 0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 2}, lie.Weight{1, 1}, lie.Weight{0, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 2}, lie.Weight{2, 1}, lie.Weight{1, 1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 1}, lie.Weight{1, 2}, lie.Weight{1, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{0, 2}, lie.Weight{0, 0}, lie.Weight{0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 1}, lie.Weight{1, 1}, lie.Weight{2, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 0, 2}, lie.Weight{0, 1, 2}, lie.Weight{0, 2, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{2, 0, 0}, lie.Weight{0, 0, 3}, lie.Weight{3, 0, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{2, 1, 0}, lie.Weight{0, 0, 1}, lie.Weight{0, 1, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 1, 1}, lie.Weight{1, 1, 1}, lie.Weight{3, 0, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 0, 1}, lie.Weight{2, 1, 0}, lie.Weight{0, 0, 3}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{1, 0, 2}, lie.Weight{1, 1, 0}, lie.Weight{1, 2, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{2, 0, 1}, lie.Weight{1, 0, 2}, lie.Weight{2, 1, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{2, 0, 1}, lie.Weight{1, 0, 2}, lie.Weight{1, 0, 1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 0, 2}, lie.Weight{1, 1, 1}, lie.Weight{1, 0, 1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 1, 1}, lie.Weight{0, 0, 1}, lie.Weight{0, 2, 0}}, 3, big.NewInt(1)},

		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{1}, lie.Weight{0}, lie.Weight{1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{2}, lie.Weight{0}, lie.Weight{3}, lie.Weight{3}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{3}, lie.Weight{3}, lie.Weight{0}, lie.Weight{0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{2}, lie.Weight{0}, lie.Weight{0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{1}, lie.Weight{2}, lie.Weight{3}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{1}, lie.Weight{2}, lie.Weight{3}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{0}, lie.Weight{1}, lie.Weight{2}, lie.Weight{0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{0}, lie.Weight{3}, lie.Weight{1}, lie.Weight{2}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{1}, lie.Weight{2}, lie.Weight{0}, lie.Weight{1}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(1), []lie.Weight{lie.Weight{2}, lie.Weight{0}, lie.Weight{0}, lie.Weight{1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 0}, lie.Weight{2, 0}, lie.Weight{2, 1}, lie.Weight{0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{0, 0}, lie.Weight{0, 1}, lie.Weight{0, 0}, lie.Weight{1, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 1}, lie.Weight{0, 2}, lie.Weight{1, 0}, lie.Weight{0, 0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 0}, lie.Weight{0, 3}, lie.Weight{0, 0}, lie.Weight{2, 0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 0}, lie.Weight{2, 1}, lie.Weight{3, 0}, lie.Weight{0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 0}, lie.Weight{2, 0}, lie.Weight{1, 0}, lie.Weight{0, 2}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{0, 2}, lie.Weight{2, 0}, lie.Weight{1, 1}, lie.Weight{0, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{1, 0}, lie.Weight{2, 1}, lie.Weight{0, 1}, lie.Weight{0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 1}, lie.Weight{2, 1}, lie.Weight{1, 0}, lie.Weight{2, 0}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(2), []lie.Weight{lie.Weight{2, 1}, lie.Weight{1, 1}, lie.Weight{2, 1}, lie.Weight{1, 0}}, 3, big.NewInt(2)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{1, 1, 0}, lie.Weight{2, 0, 1}, lie.Weight{0, 0, 2}, lie.Weight{2, 0, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 1, 1}, lie.Weight{0, 2, 1}, lie.Weight{0, 3, 0}, lie.Weight{0, 0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 2, 1}, lie.Weight{0, 1, 1}, lie.Weight{2, 0, 1}, lie.Weight{1, 1, 0}}, 3, big.NewInt(3)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 1, 1}, lie.Weight{0, 1, 0}, lie.Weight{0, 1, 1}, lie.Weight{0, 1, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{1, 0, 0}, lie.Weight{0, 0, 3}, lie.Weight{0, 0, 0}, lie.Weight{2, 0, 0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 2, 0}, lie.Weight{2, 0, 0}, lie.Weight{1, 1, 1}, lie.Weight{2, 1, 0}}, 3, big.NewInt(1)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{1, 2, 0}, lie.Weight{0, 1, 2}, lie.Weight{0, 0, 2}, lie.Weight{0, 0, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{0, 3, 0}, lie.Weight{0, 0, 3}, lie.Weight{1, 1, 1}, lie.Weight{1, 0, 2}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{1, 1, 1}, lie.Weight{0, 1, 0}, lie.Weight{2, 0, 1}, lie.Weight{1, 1, 1}}, 3, big.NewInt(0)},
		{lie.NewTypeARootSystem(3), []lie.Weight{lie.Weight{2, 0, 0}, lie.Weight{0, 0, 2}, lie.Weight{2, 0, 0}, lie.Weight{0, 0, 2}}, 3, big.NewInt(2)},
	}

	for _, c := range cases {
		alg := lie.NewAlgebra(c.rtsys)
		bun := NewCBBundle(alg, c.wts, c.ell)
		got := bun.Rank()
		if got.Cmp(c.want) != 0 {
			t.Errorf("For type A CBBundle(wts: %v, ell: %v): Rank() = %v, want %v",
				c.wts, c.ell, got, c.want)
		}
	}
}

func TestDivisor(t *testing.T) {
	cases := []struct {
		rtsys lie.RootSystem
		wt    lie.Weight
		ell   int
		n     int
		want  []*big.Rat
	}{
		// n=4, l=1
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 1, 4, []*big.Rat{big.NewRat(1, 3)}},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 1, 4, []*big.Rat{big.NewRat(2, 3)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 1, 4, []*big.Rat{big.NewRat(2, 3)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 1, 4, []*big.Rat{big.NewRat(0, 1)}},

		// n=4, l=2
		{lie.NewTypeARootSystem(1), lie.Weight{0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(1), lie.Weight{1}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(1), lie.Weight{2}, 2, 4, []*big.Rat{big.NewRat(2, 3)}},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 2, 4, []*big.Rat{big.NewRat(1, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 1}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 1}, 2, 4, []*big.Rat{big.NewRat(2, 3)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 2, 4, []*big.Rat{big.NewRat(4, 3)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 2, 4, []*big.Rat{big.NewRat(1, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 0}, 2, 4, []*big.Rat{big.NewRat(2, 3)}},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},

		// n=6
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 2, 6, []*big.Rat{big.NewRat(0, 1), big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 2, 6, []*big.Rat{big.NewRat(1, 1), big.NewRat(3, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 2, 6, []*big.Rat{big.NewRat(8, 5), big.NewRat(4, 5)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 2, 6, []*big.Rat{big.NewRat(27, 5), big.NewRat(31, 5)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 0}, 3, 6, []*big.Rat{big.NewRat(0, 1), big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 0, 2}, 3, 6, []*big.Rat{big.NewRat(0, 1), big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 0}, 3, 6, []*big.Rat{big.NewRat(0, 1), big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 1, 2}, 3, 6, []*big.Rat{big.NewRat(4, 1), big.NewRat(8, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 2, 0}, 3, 6, []*big.Rat{big.NewRat(64, 5), big.NewRat(62, 5)}},
		{lie.NewTypeARootSystem(3), lie.Weight{0, 3, 0}, 3, 6, []*big.Rat{big.NewRat(12, 5), big.NewRat(6, 5)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 0, 1}, 3, 6, []*big.Rat{big.NewRat(22, 1), big.NewRat(34, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{1, 1, 1}, 3, 6, []*big.Rat{big.NewRat(122, 1), big.NewRat(137, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 0, 0}, 3, 6, []*big.Rat{big.NewRat(0, 1), big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(3), lie.Weight{2, 1, 0}, 3, 6, []*big.Rat{big.NewRat(4, 1), big.NewRat(8, 1)}},

		// Allowing zero rank
		{lie.NewTypeARootSystem(2), lie.Weight{0, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 1}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{0, 2}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{1, 1}, 2, 4, []*big.Rat{big.NewRat(1, 1)}},
		{lie.NewTypeARootSystem(2), lie.Weight{2, 0}, 2, 4, []*big.Rat{big.NewRat(0, 1)}},
	}

	for _, c := range cases {
		alg := lie.NewAlgebra(c.rtsys)
		bun := NewSymmetricCBBundle(alg, c.wt, c.ell, c.n)
		got := bun.SymmetrizedDivisor()
		if len(got) != len(c.want) {
			t.Errorf("For type A CBBundle(wt: %v, ell: %v, n: %v): wrong divisor dimension %v, want %v",
				c.wt, c.ell, c.n, len(got), len(c.want))
		}
		for i := range got {
			if got[i].Cmp(c.want[i]) != 0 {
				t.Errorf("For type A CBBundle(wt: %v, ell: %v, n: %v): Divisor() = %v, want %v",
					c.wt, c.ell, c.n, got, c.want)
			}
		}
	}
}

func BenchmarkSymmetricCBRank(b *testing.B) {
	rank := 5
	level := 4
	n := 10
	alg := lie.NewAlgebra(lie.NewTypeARootSystem(rank))
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < 1; i++ {
		for j := 0; j < len(wts); j++ {
			bun := NewSymmetricCBBundle(alg, wts[j], level, n)
			bun.Rank()
			// if rk.Cmp(big.NewInt(0)) == 0 {
			// 	continue
			// }
			// fmt.Printf("%v: %v\n", wts[j], rk)
		}
	}
}

func BenchmarkSymmetricDivisor(b *testing.B) {
	rank := 5
	level := 4
	n := 100
	alg := lie.NewAlgebra(lie.NewTypeARootSystem(rank))
	wts := alg.Weights(level)

	b.ReportAllocs()
	for i := 0; i < 1; i++ {
		for j := 0; j < len(wts); j++ {
			bun := NewSymmetricCBBundle(alg, wts[j], level, n)
			// rk := bun.Rank()
			bun.SymmetrizedDivisor()
			// if rk.Cmp(big.NewInt(0)) == 0 {
			// 	continue
			// }
			// fmt.Printf("%v: %v %v\n", wts[j], rk, div)
		}
	}
}
