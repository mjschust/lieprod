package lie

import (
	"testing"

	"github.com/mjschust/cblocks/util"
)

func TestTypeAConvertWeight(t *testing.T) {
	cases := []struct {
		alg  TypeA
		wt   Weight
		want epCoord
	}{
		{TypeA{1}, Weight{0}, epCoord{0, 0}},
		{TypeA{1}, Weight{1}, epCoord{1, 0}},
		{TypeA{1}, Weight{2}, epCoord{2, 0}},
		{TypeA{2}, Weight{0, 0}, epCoord{0, 0, 0}},
		{TypeA{2}, Weight{1, 0}, epCoord{1, 0, 0}},
		{TypeA{2}, Weight{0, 1}, epCoord{1, 1, 0}},
		{TypeA{2}, Weight{1, 1}, epCoord{2, 1, 0}},
	}

	for _, c := range cases {
		got := c.alg.convertWeight(c.wt)
		if !equals(got, c.want) {
			t.Errorf("convertWeight(%v) = %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeAConvertEpCoord(t *testing.T) {
	cases := []struct {
		alg  TypeA
		epc  epCoord
		want Weight
	}{
		{TypeA{1}, epCoord{0, 0}, Weight{0}},
		{TypeA{1}, epCoord{1, 1}, Weight{0}},
		{TypeA{1}, epCoord{1, 0}, Weight{1}},
		{TypeA{1}, epCoord{0, 1}, Weight{-1}},
		{TypeA{2}, epCoord{0, 0, 0}, Weight{0, 0}},
		{TypeA{2}, epCoord{1, 1, 1}, Weight{0, 0}},
		{TypeA{2}, epCoord{1, 0, 0}, Weight{1, 0}},
		{TypeA{2}, epCoord{1, 1, 0}, Weight{0, 1}},
		{TypeA{2}, epCoord{2, 1, 0}, Weight{1, 1}},
		{TypeA{2}, epCoord{1, 2, 0}, Weight{-1, 2}},
	}

	for _, c := range cases {
		got := c.alg.convertEpCoord(c.epc)
		if !equals(got, c.want) {
			t.Errorf("convertEpCoord(%v) = %v, want %v", c.epc, got, c.want)
		}
	}
}

func TestTypeADualCoxeter(t *testing.T) {
	cases := []struct {
		alg  TypeA
		want int
	}{
		{TypeA{1}, 2},
		{TypeA{2}, 3},
		{TypeA{3}, 4},
	}

	for _, c := range cases {
		got := c.alg.DualCoxeter()
		if got != c.want {
			t.Errorf("DualCoxeter() == %v, want %v", got, c.want)
		}
	}
}

func TestTypeAPositiveRoots(t *testing.T) {
	cases := []struct {
		alg  TypeA
		want []Root
	}{
		{TypeA{1}, []Root{Root{1}}},
		{TypeA{2}, []Root{Root{1, 0}, Root{1, 1}, Root{0, 1}}},
		{TypeA{3}, []Root{
			Root{1, 0, 0},
			Root{1, 1, 0},
			Root{1, 1, 1},
			Root{0, 1, 0},
			Root{0, 1, 1},
			Root{0, 0, 1},
		}},
	}

	for _, c := range cases {
		got := c.alg.PositiveRoots()
		if len(got) != len(c.want) {
			t.Errorf("PositiveRoots() == %v, want %v", got, c.want)
		}
		for i := range c.want {
			if !equals(got[i], c.want[i]) {
				t.Errorf("PositiveRoots() == %v, want %v", got, c.want)
			}
		}
	}
}

func TestTypeAKillingForm(t *testing.T) {
	cases := []struct {
		alg      TypeA
		wt1, wt2 Weight
		want     float64
	}{
		{TypeA{1}, Weight{0}, Weight{0}, 0},
		{TypeA{1}, Weight{1}, Weight{0}, 0},
		{TypeA{1}, Weight{0}, Weight{1}, 0},
		{TypeA{1}, Weight{1}, Weight{1}, 0.5},
		{TypeA{1}, Weight{2}, Weight{1}, 1},
		{TypeA{1}, Weight{1}, Weight{2}, 1},
		{TypeA{1}, Weight{2}, Weight{2}, 2},
		{TypeA{2}, Weight{0, 0}, Weight{0, 0}, 0},
		{TypeA{2}, Weight{1, 0}, Weight{0, 0}, 0},
		{TypeA{2}, Weight{0, 0}, Weight{1, 0}, 0},
		{TypeA{2}, Weight{1, 0}, Weight{1, 0}, 0.6666666666666667},
		{TypeA{2}, Weight{0, 1}, Weight{1, 0}, 0.33333333333333337},
		{TypeA{2}, Weight{1, 0}, Weight{0, 1}, 0.33333333333333337},
		{TypeA{2}, Weight{0, 1}, Weight{0, 1}, 0.6666666666666667},
	}

	for _, c := range cases {
		got := c.alg.KillingForm(c.wt1, c.wt2)
		if got != c.want {
			t.Errorf("KillingForm(%v, %v) == %v, want %v", c.wt1, c.wt2, got, c.want)
		}
	}
}

func TestTypeALevel(t *testing.T) {
	cases := []struct {
		alg  TypeA
		wt   Weight
		want float64
	}{
		{TypeA{1}, Weight{0}, 0},
		{TypeA{1}, Weight{1}, 1},
		{TypeA{1}, Weight{2}, 2},
		{TypeA{2}, Weight{0, 0}, 0},
		{TypeA{2}, Weight{1, 0}, 1},
		{TypeA{2}, Weight{0, 1}, 1},
		{TypeA{2}, Weight{1, 1}, 2},
	}

	for _, c := range cases {
		got := c.alg.Level(c.wt)
		if got != c.want {
			t.Errorf("Level(%v) == %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeADual(t *testing.T) {
	cases := []struct {
		alg      TypeA
		wt, want Weight
	}{
		{TypeA{1}, Weight{0}, Weight{0}},
		{TypeA{1}, Weight{1}, Weight{1}},
		{TypeA{1}, Weight{2}, Weight{2}},
		{TypeA{2}, Weight{0, 0}, Weight{0, 0}},
		{TypeA{2}, Weight{1, 0}, Weight{0, 1}},
		{TypeA{2}, Weight{0, 1}, Weight{1, 0}},
		{TypeA{2}, Weight{1, 1}, Weight{1, 1}},
	}

	for _, c := range cases {
		got := c.alg.Dual(c.wt)
		if !equals(got, c.want) {
			t.Errorf("Dual(%v) = %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeAReflectIntoChamber(t *testing.T) {
	cases := []struct {
		alg      TypeA
		wt, want Weight
		parity   int
	}{
		{TypeA{1}, Weight{0}, Weight{0}, 1},
		{TypeA{1}, Weight{1}, Weight{1}, 1},
		{TypeA{1}, Weight{-1}, Weight{1}, -1},
		{TypeA{2}, Weight{0, 0}, Weight{0, 0}, 1},
		{TypeA{2}, Weight{1, 0}, Weight{1, 0}, 1},
		{TypeA{2}, Weight{0, 1}, Weight{0, 1}, 1},
		{TypeA{2}, Weight{-1, 0}, Weight{0, 1}, 1},
		{TypeA{2}, Weight{0, -1}, Weight{1, 0}, 1},
		{TypeA{2}, Weight{-1, -1}, Weight{1, 1}, -1},
	}

	for _, c := range cases {
		got, parity := c.alg.ReflectToChamber(c.wt)
		if !equals(got, c.want) || parity != c.parity {
			t.Errorf("ReflectToChamber(%v) = %v, %v, want %v, %v",
				c.wt, got, parity, c.want, c.parity)
		}
	}
}

func TestTypeAOrbitIterator(t *testing.T) {
	cases := []struct {
		alg   TypeA
		wt    Weight
		orbit []Weight
	}{
		{TypeA{1}, Weight{0}, []Weight{Weight{0}}},
		{TypeA{1}, Weight{1}, []Weight{Weight{1}, Weight{-1}}},
		{TypeA{1}, Weight{2}, []Weight{Weight{2}, Weight{-2}}},
		{TypeA{2}, Weight{0, 0}, []Weight{Weight{0, 0}}},
		{TypeA{2}, Weight{1, 0}, []Weight{
			Weight{1, 0},
			Weight{-1, 1},
			Weight{0, -1}}},
		{TypeA{2}, Weight{0, 1}, []Weight{
			Weight{0, 1},
			Weight{1, -1},
			Weight{-1, 0}}},
		{TypeA{2}, Weight{1, 1}, []Weight{
			Weight{1, 1},
			Weight{-1, 2},
			Weight{2, -1},
			Weight{1, -2},
			Weight{-2, 1},
			Weight{-1, -1}}},
	}

	for _, c := range cases {
		orbitSet := weightSetFromList(c.orbit)
		orbitIter := c.alg.NewOrbitIterator(c.wt)
		orbitSize := 0
		for orbitIter.HasNext() {
			nextWt := orbitIter.Next()
			_, present := orbitSet.Get(nextWt)
			if !present {
				t.Errorf("OrbitIterator(%v) does not contain %v", c.wt, nextWt)
			}
			orbitSize++
		}
		if orbitSize != len(c.orbit) {
			t.Errorf("OrbitIterator(%v) is missing orbit elements", c.wt)
		}
	}
}

func weightSetFromList(wts []Weight) util.VectorMap {
	vmap := util.NewVectorMap()
	for _, wt := range wts {
		vmap.Put(wt, true)
	}

	return vmap
}

func equals(v1, v2 []float64) bool {
	if (v1 == nil) != (v2 == nil) {
		return false
	}

	if len(v1) != len(v2) {
		return false
	}

	for i := range v1 {
		if v1[i] != v2[i] {
			return false
		}
	}

	return true
}
