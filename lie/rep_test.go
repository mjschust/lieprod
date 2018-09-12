package lie

import "testing"

func TestReprDimension(t *testing.T) {
	cases := []struct {
		alg       TypeA
		highestWt Weight
		want      float64
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
