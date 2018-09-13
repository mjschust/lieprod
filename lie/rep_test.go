package lie

import (
	"testing"
)

func TestReprDimension(t *testing.T) {
	cases := []struct {
		alg       TypeA
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
		got := ReprDimension(c.alg, c.highestWt)
		if got != c.want {
			t.Errorf("ReprDimension(%v, %v) = %v, want %v", c.alg, c.highestWt, got, c.want)
		}
	}
}

func TestDominantChar(t *testing.T) {
	cases := []struct {
		alg       TypeA
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
		domChar := DominantChar(c.alg, c.highestWt)
		for i := range c.wantWts {
			gotMult, present := domChar.Get(c.wantWts[i])
			if !present {
				t.Errorf("DominantChar(%v) missing weight %v", c.highestWt, c.wantWts[i])
				continue
			}
			if gotMult.(int) != c.wantMults[i] {
				t.Errorf("DominantChar(%v)[%v] = %v, want %v", c.highestWt, c.wantWts[i], gotMult, c.wantMults[i])
			}
		}
	}
}

func TestTensor(t *testing.T) {
	cases := []struct {
		alg       TypeA
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
		tensorDecomp := Tensor(c.alg, c.wt1, c.wt2)
		for i := range c.wantWts {
			gotMult, present := tensorDecomp.Get(c.wantWts[i])
			if !present {
				t.Errorf("Tensor(%v, %v) missing weight %v", c.wt1, c.wt2, c.wantWts[i])
				continue
			}
			if gotMult.(int) != c.wantMults[i] {
				t.Errorf("Tensor(%v, %v)[%v] = %v, want %v", c.wt1, c.wt2, c.wantWts[i], gotMult, c.wantMults[i])
			}
			tensorDecomp.Remove(c.wantWts[i])
		}
		if tensorDecomp.Size() != 0 {
			t.Errorf("Tensor(%v, %v) contains extra weights", c.wt1, c.wt2)
		}
	}
}

func BenchmarkTensor(b *testing.B) {
	rank := 4
	level := 4
	alg := TypeA{rank}
	wts := alg.Weights(level)
	for i := 0; i < b.N; i++ {
		for i := range wts {
			for j := range wts {
				Tensor(alg, wts[i], wts[j])
			}
		}
	}
}
