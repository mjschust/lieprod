package bundle

import (
	"testing"

	"github.com/mjschust/cblocks/lie"
)

func TestSymmetricRank(t *testing.T) {

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
