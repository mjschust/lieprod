package lie

import (
	"math/big"
	"math/rand"
	"testing"
)

func TestReprDimension(t *testing.T) {
	cases := []struct {
		rtsys     RootSystem
		highestWt Weight
		want      int
	}{
		{TypeA{1}, Weight{0}, 1},
		{TypeA{1}, Weight{1}, 2},
		{TypeA{1}, Weight{2}, 3},
		{TypeA{2}, Weight{0, 0}, 1},
		{TypeA{2}, Weight{1, 0}, 3},
		{TypeA{2}, Weight{0, 1}, 3},
		{TypeA{2}, Weight{1, 1}, 8},
		{TypeA{2}, Weight{2, 1}, 15},
		{TypeA{3}, Weight{0, 0, 0}, 1},
		{TypeA{3}, Weight{1, 0, 0}, 4},
		{TypeA{3}, Weight{0, 1, 0}, 6},
		{TypeA{3}, Weight{0, 0, 1}, 4},
	}

	for _, c := range cases {
		alg := NewAlgebra(c.rtsys)
		got := alg.ReprDimension(c.highestWt)
		if got.Cmp(big.NewInt(int64(c.want))) != 0 {
			t.Errorf("ReprDimension(%v, %v) = %v, want %v", c.rtsys, c.highestWt, got, c.want)
		}
	}
}

func TestDominantChar(t *testing.T) {
	cases := []struct {
		rtsys     RootSystem
		highestWt Weight
		wantWts   [][]int
		wantMults []int
	}{
		{TypeA{1}, Weight{0}, [][]int{{0}}, []int{1}},
		{TypeA{1}, Weight{1}, [][]int{{1}}, []int{1}},
		{TypeA{1}, Weight{2}, [][]int{{0}, {2}}, []int{1, 1}},
		{TypeA{1}, Weight{3}, [][]int{{1}, {3}}, []int{1, 1}},
		{TypeA{1}, Weight{4}, [][]int{{0}, {2}, {4}}, []int{1, 1, 1}},
		{TypeA{2}, Weight{0, 0}, [][]int{{0, 0}}, []int{1}},
		{TypeA{2}, Weight{1, 0}, [][]int{{1, 0}}, []int{1}},
		{TypeA{2}, Weight{0, 1}, [][]int{{0, 1}}, []int{1}},
		{TypeA{2}, Weight{1, 1}, [][]int{{1, 1}, {0, 0}}, []int{1, 2}},
		{TypeA{2}, Weight{2, 1}, [][]int{{2, 1}, {1, 0}, {0, 2}}, []int{1, 2, 1}},
		{TypeA{2}, Weight{2, 3},
			[][]int{
				{0, 1},
				{0, 4},
				{1, 2},
				{2, 0},
				{2, 3},
				{3, 1}},
			[]int{3, 1, 2, 2, 1, 1}},
		{TypeA{3}, Weight{1, 2, 1},
			[][]int{
				{0, 0, 0},
				{0, 1, 2},
				{0, 2, 0},
				{1, 0, 1},
				{1, 2, 1},
				{2, 0, 2},
				{2, 1, 0}},
			[]int{7, 2, 4, 5, 1, 1, 2}},
	}

	for _, c := range cases {
		alg := NewAlgebra(c.rtsys)
		domChar := alg.DominantChar(c.highestWt)
		for i := range c.wantWts {
			gotMult := domChar.Multiplicity(c.wantWts[i])
			if gotMult.Cmp(big.NewInt(int64(c.wantMults[i]))) != 0 {
				t.Errorf("DominantChar(%v)[%v] = %v, want %v", c.highestWt, c.wantWts[i], gotMult, c.wantMults[i])
			}
		}
	}
}

func TestTensor(t *testing.T) {
	cases := []struct {
		rtsys     RootSystem
		wt1, wt2  Weight
		wantWts   [][]int
		wantMults []int
	}{
		{TypeA{1}, Weight{0}, Weight{0},
			[][]int{{0}},
			[]int{1}},
		{TypeA{1}, Weight{1}, Weight{0},
			[][]int{{1}},
			[]int{1}},
		{TypeA{1}, Weight{0}, Weight{1},
			[][]int{{1}},
			[]int{1}},
		{TypeA{1}, Weight{1}, Weight{1},
			[][]int{{2}, {0}},
			[]int{1, 1}},
		{TypeA{1}, Weight{2}, Weight{1},
			[][]int{{3}, {1}},
			[]int{1, 1}},
		{TypeA{1}, Weight{2}, Weight{2},
			[][]int{{4}, {2}, {0}},
			[]int{1, 1, 1}},
		{TypeA{2}, Weight{0, 0}, Weight{0, 0},
			[][]int{{0, 0}},
			[]int{1}},
		{TypeA{2}, Weight{1, 0}, Weight{0, 0},
			[][]int{{1, 0}},
			[]int{1}},
		{TypeA{2}, Weight{0, 0}, Weight{0, 1},
			[][]int{{0, 1}},
			[]int{1}},
		{TypeA{2}, Weight{1, 0}, Weight{1, 0},
			[][]int{{0, 1}, {2, 0}},
			[]int{1, 1}},
		{TypeA{2}, Weight{1, 0}, Weight{0, 1},
			[][]int{{0, 0}, {1, 1}},
			[]int{1, 1}},
		{TypeA{2}, Weight{1, 1}, Weight{0, 1},
			[][]int{{0, 1}, {1, 2}, {2, 0}},
			[]int{1, 1, 1}},
		{TypeA{2}, Weight{1, 1}, Weight{1, 1},
			[][]int{{0, 0}, {0, 3}, {3, 0}, {1, 1}, {2, 2}},
			[]int{1, 1, 1, 2, 1}},
		{TypeA{2}, Weight{2, 1}, Weight{1, 1},
			[][]int{{0, 2}, {1, 0}, {1, 3}, {2, 1}, {3, 2}, {4, 0}},
			[]int{1, 1, 1, 2, 1, 1}},
		{TypeA{3}, Weight{1, 0, 1}, Weight{0, 2, 1},
			[][]int{{0, 1, 3}, {0, 2, 1}, {1, 0, 2}, {1, 1, 0}, {1, 2, 2}, {1, 3, 0}, {2, 1, 1}},
			[]int{1, 2, 1, 1, 1, 1, 1}},
	}

	for _, c := range cases {
		alg := algebraImpl{c.rtsys}
		tensorDecomp := alg.tensorProduct(c.wt1, c.wt2)
		if len(tensorDecomp.Weights()) != len(c.wantWts) {
			t.Errorf("Tensor(%v, %v) contains wrong number of weights", c.wt1, c.wt2)
		}
		for i := range c.wantWts {
			gotMult := tensorDecomp.Multiplicity(c.wantWts[i])
			if gotMult.Cmp(big.NewInt(int64(c.wantMults[i]))) != 0 {
				t.Errorf("Tensor(%v, %v)[%v] = %v, want %v", c.wt1, c.wt2, c.wantWts[i], gotMult, c.wantMults[i])
			}
		}
	}
}

func TestFusion(t *testing.T) {
	cases := []struct {
		rtsys     RootSystem
		ell       int
		wt1, wt2  Weight
		wantWts   [][]int
		wantMults []int
	}{
		{TypeA{1}, 1, Weight{0}, Weight{0},
			[][]int{{0}},
			[]int{1}},
		{TypeA{1}, 1, Weight{1}, Weight{0},
			[][]int{{1}},
			[]int{1}},
		{TypeA{1}, 1, Weight{0}, Weight{1},
			[][]int{{1}},
			[]int{1}},
		{TypeA{1}, 1, Weight{1}, Weight{1},
			[][]int{{0}},
			[]int{1}},
		{TypeA{1}, 2, Weight{1}, Weight{1},
			[][]int{{0}, {2}},
			[]int{1, 1}},
		{TypeA{1}, 2, Weight{2}, Weight{1},
			[][]int{{1}},
			[]int{1}},
		{TypeA{1}, 2, Weight{2}, Weight{2},
			[][]int{{0}},
			[]int{1}},
		{TypeA{2}, 1, Weight{0, 1}, Weight{0, 1},
			[][]int{{1, 0}},
			[]int{1}},
		{TypeA{2}, 1, Weight{0, 1}, Weight{1, 0},
			[][]int{{0, 0}},
			[]int{1}},
		{TypeA{2}, 1, Weight{1, 0}, Weight{1, 0},
			[][]int{{0, 1}},
			[]int{1}},
		{TypeA{2}, 2, Weight{1, 1}, Weight{0, 1},
			[][]int{{2, 0}, {0, 1}},
			[]int{1, 1}},
		{TypeA{2}, 4, Weight{1, 1}, Weight{1, 1},
			[][]int{{2, 2}, {3, 0}, {0, 3}, {0, 0}, {1, 1}},
			[]int{1, 1, 1, 1, 2}},
		{TypeA{2}, 4, Weight{2, 2}, Weight{2, 2},
			[][]int{{1, 1}, {0, 0}, {2, 2}},
			[]int{1, 1, 1}},
	}

	for _, c := range cases {
		alg := algebraImpl{c.rtsys}
		fusionDecomp := alg.fusionProduct(c.ell, c.wt1, c.wt2)
		// if len(fusionDecomp.Weights()) != len(c.wantWts) {
		// 	t.Errorf("Fusion(%v, %v, %v) contains wrong number of weights", c.ell, c.wt1, c.wt2)
		// }
		for i := range c.wantWts {
			gotMult := fusionDecomp.Multiplicity(c.wantWts[i])
			if gotMult.Cmp(big.NewInt(int64(c.wantMults[i]))) != 0 {
				t.Errorf("Fusion(%v, %v, %v)[%v] = %v, want %v", c.ell, c.wt1, c.wt2, c.wantWts[i], gotMult, c.wantMults[i])
			}
		}
	}
}

func TestMultiTensor(t *testing.T) {
	cases := []struct {
		rtsys     RootSystem
		wts       []Weight
		wantWts   [][]int
		wantMults []int
	}{
		{
			TypeA{1},
			[]Weight{Weight{0}, Weight{0}, Weight{0}},
			[][]int{{0}},
			[]int{1},
		},
		{
			TypeA{1},
			[]Weight{Weight{1}, Weight{0}, Weight{0}},
			[][]int{{1}},
			[]int{1},
		},
		{
			TypeA{1},
			[]Weight{Weight{1}, Weight{0}, Weight{1}},
			[][]int{{2}, {0}},
			[]int{1, 1},
		},
		{
			TypeA{1},
			[]Weight{Weight{1}, Weight{1}, Weight{1}},
			[][]int{{3}, {1}},
			[]int{1, 2},
		},
		{
			TypeA{1},
			[]Weight{Weight{1}, Weight{1}, Weight{2}},
			[][]int{{4}, {2}, {0}},
			[]int{1, 2, 1},
		},
		{
			TypeA{2},
			[]Weight{Weight{0, 0}, Weight{0, 0}, Weight{0, 0}},
			[][]int{{0, 0}},
			[]int{1},
		},
		{
			TypeA{2},
			[]Weight{Weight{0, 0}, Weight{2, 0}, Weight{0, 0}},
			[][]int{{2, 0}},
			[]int{1},
		},
		{
			TypeA{2},
			[]Weight{Weight{0, 0}, Weight{1, 0}, Weight{0, 1}},
			[][]int{{0, 0}, {1, 1}},
			[]int{1, 1},
		},
		{
			TypeA{2},
			[]Weight{Weight{1, 0}, Weight{1, 1}, Weight{0, 1}},
			[][]int{{0, 0}, {0, 3}, {1, 1}, {2, 2}, {3, 0}},
			[]int{1, 1, 3, 1, 1},
		},
		{
			TypeA{2},
			[]Weight{Weight{1, 0}, Weight{1, 1}, Weight{2, 0}},
			[][]int{{0, 0}, {0, 3}, {1, 1}, {2, 2}, {3, 0}, {4, 1}},
			[]int{1, 1, 3, 2, 2, 1},
		},
		{
			TypeA{3},
			[]Weight{Weight{1, 0, 1}, Weight{0, 2, 1}, Weight{1, 1, 0}, Weight{2, 0, 0}},
			[][]int{{0, 0, 0}, {0, 0, 4}, {0, 1, 2}, {0, 2, 0}, {0, 2, 4}, {0, 3, 2},
				{0, 4, 0}, {0, 5, 2}, {0, 6, 0}, {1, 0, 1}, {1, 0, 5}, {1, 1, 3},
				{1, 2, 1}, {1, 3, 3}, {1, 4, 1}, {2, 0, 2}, {2, 1, 0}, {2, 1, 4},
				{2, 2, 2}, {2, 3, 0}, {2, 4, 2}, {2, 5, 0}, {3, 0, 3}, {3, 1, 1},
				{3, 2, 3}, {3, 3, 1}, {4, 0, 0}, {4, 0, 4}, {4, 1, 2}, {4, 2, 0},
				{4, 3, 2}, {4, 4, 0}, {5, 0, 1}, {5, 1, 3}, {5, 2, 1}, {6, 0, 2},
				{6, 1, 0}},
			[]int{5, 10, 36, 24, 6, 21, 15, 1, 1, 29, 5, 38, 57, 6, 14, 43, 33, 9,
				34, 25, 2, 2, 22, 39, 5, 12, 8, 3, 14, 11, 1, 1, 6, 1, 3, 1, 1},
		},
	}

	for _, c := range cases {
		alg := NewAlgebra(c.rtsys)
		tensorDecomp := alg.Tensor(c.wts...)
		if len(tensorDecomp.Weights()) != len(c.wantWts) {
			t.Errorf("Tensor(%v) contains wrong number of weights", c.wts)
		}
		for i := range c.wantWts {
			gotMult := tensorDecomp.Multiplicity(c.wantWts[i])
			if gotMult.Cmp(big.NewInt(int64(c.wantMults[i]))) != 0 {
				t.Errorf("Tensor(%v)[%v] = %v, want %v", c.wts, c.wantWts[i], gotMult, c.wantMults[i])
			}
		}

	}
}

func BenchmarkTensorSmall(b *testing.B) {
	rank := 4
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		for j := range wts {
			for k := range wts {
				alg.tensorProduct(wts[j], wts[k])
			}
		}
	}
}

func BenchmarkMultiTensorSmall(b *testing.B) {
	rank := 3
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		for j := range wts {
			alg.Tensor(wts[j], wts[j], wts[j])
		}

	}
}

func BenchmarkMultiTensorLarge(b *testing.B) {
	rank := 5
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		j := rand.Intn(len(wts))
		k := rand.Intn(len(wts))
		l := rand.Intn(len(wts))
		alg.Tensor(wts[j], wts[k], wts[l])
	}
}

func BenchmarkTensorLarge(b *testing.B) {
	rank := 6
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		j := rand.Intn(len(wts))
		k := rand.Intn(len(wts))
		alg.tensorProduct(wts[j], wts[k])
	}
}

func BenchmarkTensorParallel(b *testing.B) {
	numRoutines := 100
	rank := 6
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		done := make(chan int, numRoutines)
		for j := 0; j < numRoutines; j++ {
			go func(c chan int) {
				alg.tensorProduct(wts[rand.Intn(len(wts))], wts[rand.Intn(len(wts))])
				c <- 0
			}(done)
		}

		for j := 0; j < numRoutines; j++ {
			<-done
		}
	}
}

func BenchmarkCBRank(b *testing.B) {
	rank := 5
	level := 4
	alg := algebraImpl{TypeA{rank}}
	wts := alg.Weights(level)
	b.ReportAllocs()
	for i := 0; i < 1; i++ {
		for j := 0; j < len(wts); j++ {
			wt10 := []Weight{wts[j], wts[j], wts[j], wts[j], wts[j], wts[j], wts[j], wts[j], wts[j], wts[j]}
			alg.CBRank(level, wt10...)
			// if rk.Cmp(big.NewInt(0)) == 0 {
			// 	continue
			// }
			// fmt.Printf("%v: %v\n", wts[j], rk)
			// prod := alg.fusionProduct(level, wts[j], wts[j])
			// fmt.Printf("%v: ", wts[j])
			// for _, wt := range prod.Weights() {
			// 	fmt.Printf("(%v, %v), ", wt, prod.Multiplicity(wt))
			// }
			// fmt.Printf("\n")
		}
	}
}
