package bundle

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