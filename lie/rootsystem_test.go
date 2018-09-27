package lie

import (
	"testing"

	"github.com/mjschust/cblocks/util"
)

func TestTypeAConvertWeightToEpc(t *testing.T) {
	cases := []struct {
		rtsys typeA
		wt    Weight
		want  []int
	}{
		{typeA{1}, Weight{0}, []int{0, 0}},
		{typeA{1}, Weight{1}, []int{1, 0}},
		{typeA{1}, Weight{2}, []int{2, 0}},
		{typeA{2}, Weight{0, 0}, []int{0, 0, 0}},
		{typeA{2}, Weight{1, 0}, []int{1, 0, 0}},
		{typeA{2}, Weight{0, 1}, []int{1, 1, 0}},
		{typeA{2}, Weight{1, 1}, []int{2, 1, 0}},
	}

	for _, c := range cases {
		got := make([]int, c.rtsys.rank+1)
		c.rtsys.convertWeightToEpc(c.wt, got)
		if !equals(got, c.want) {
			t.Errorf("convertWeightToEpc(%v) = %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeAConvertEpCoord(t *testing.T) {
	cases := []struct {
		rtsys typeA
		epc   []int
		want  Weight
	}{
		{typeA{1}, []int{0, 0}, Weight{0}},
		{typeA{1}, []int{1, 1}, Weight{0}},
		{typeA{1}, []int{1, 0}, Weight{1}},
		{typeA{1}, []int{0, 1}, Weight{-1}},
		{typeA{2}, []int{0, 0, 0}, Weight{0, 0}},
		{typeA{2}, []int{1, 1, 1}, Weight{0, 0}},
		{typeA{2}, []int{1, 0, 0}, Weight{1, 0}},
		{typeA{2}, []int{1, 1, 0}, Weight{0, 1}},
		{typeA{2}, []int{2, 1, 0}, Weight{1, 1}},
		{typeA{2}, []int{1, 2, 0}, Weight{-1, 2}},
	}

	for _, c := range cases {
		var got Weight = make([]int, c.rtsys.rank)
		c.rtsys.convertEpCoord(c.epc, got)
		if !equals(got, c.want) {
			t.Errorf("convertEpCoord(%v) = %v, want %v", c.epc, got, c.want)
		}
	}
}

func TestTypeADualCoxeter(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		want  int
	}{
		{typeA{1}, 2},
		{typeA{2}, 3},
		{typeA{3}, 4},
	}

	for _, c := range cases {
		got := c.rtsys.DualCoxeter()
		if got != c.want {
			t.Errorf("DualCoxeter() == %v, want %v", got, c.want)
		}
	}
}

func TestTypeAPositiveRoots(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		want  []Root
	}{
		{typeA{1}, []Root{Root{1}}},
		{typeA{2}, []Root{Root{1, 0}, Root{1, 1}, Root{0, 1}}},
		{typeA{3}, []Root{
			Root{1, 0, 0},
			Root{1, 1, 0},
			Root{1, 1, 1},
			Root{0, 1, 0},
			Root{0, 1, 1},
			Root{0, 0, 1},
		}},
	}

	for _, c := range cases {
		got := c.rtsys.PositiveRoots()
		if len(got) != len(c.want) {
			t.Errorf("len(PositiveRoots()) == %v, want %v", len(got), len(c.want))
		}
		for i := range c.want {
			if !equals(got[i], c.want[i]) {
				t.Errorf("PositiveRoots() == %v, want %v", got, c.want)
			}
		}
	}
}

func TestTypeAWeights(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		level int
		want  [][]int
	}{
		{typeA{1}, 0, [][]int{{0}}},
		{typeA{1}, 1, [][]int{{0}, {1}}},
		{typeA{1}, 2, [][]int{{0}, {1}, {2}}},
		{typeA{2}, 0, [][]int{{0, 0}}},
		{typeA{2}, 1, [][]int{{0, 0}, {1, 0}, {0, 1}}},
		{typeA{2}, 2, [][]int{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 0}, {0, 2}}},
		{typeA{3}, 0, [][]int{
			{0, 0, 0},
		}},
		{typeA{3}, 1, [][]int{
			{0, 0, 0},
			{1, 0, 0},
			{0, 1, 0},
			{0, 0, 1},
		}},
		{typeA{3}, 2, [][]int{
			{0, 0, 0},
			{1, 0, 0},
			{0, 1, 0},
			{0, 0, 1},
			{1, 1, 0},
			{0, 1, 1},
			{1, 0, 1},
			{2, 0, 0},
			{0, 2, 0},
			{0, 0, 2},
		}},
	}

	for _, c := range cases {
		wantSet := util.NewVectorMap()
		for _, wt := range c.want {
			wantSet.Put(wt, true)
		}
		got := c.rtsys.Weights(c.level)
		if len(got) != len(c.want) {
			t.Errorf("len(Weights(%v)) == %v, want %v", c.level, len(got), len(c.want))
		}
		for _, gotWt := range got {
			_, present := wantSet.Remove(gotWt)
			if !present {
				t.Errorf("Weights(%v) should not contain %v", c.level, gotWt)
			}
		}
		if wantSet.Size() != 0 {
			t.Errorf("Weight(%v) is missing %v", c.level, wantSet.Keys())
		}
	}
}

func TestTypeAKillingForm(t *testing.T) {
	cases := []struct {
		rtsys    RootSystem
		wt1, wt2 Weight
		want     float64
	}{
		{typeA{1}, Weight{0}, Weight{0}, 0},
		{typeA{1}, Weight{1}, Weight{0}, 0},
		{typeA{1}, Weight{0}, Weight{1}, 0},
		{typeA{1}, Weight{1}, Weight{1}, 0.5},
		{typeA{1}, Weight{2}, Weight{1}, 1},
		{typeA{1}, Weight{1}, Weight{2}, 1},
		{typeA{1}, Weight{2}, Weight{2}, 2},
		{typeA{2}, Weight{0, 0}, Weight{0, 0}, 0},
		{typeA{2}, Weight{1, 0}, Weight{0, 0}, 0},
		{typeA{2}, Weight{0, 0}, Weight{1, 0}, 0},
		{typeA{2}, Weight{1, 0}, Weight{1, 0}, 0.6666666666666666},
		{typeA{2}, Weight{0, 1}, Weight{1, 0}, 0.33333333333333333},
		{typeA{2}, Weight{1, 0}, Weight{0, 1}, 0.33333333333333333},
		{typeA{2}, Weight{0, 1}, Weight{0, 1}, 0.6666666666666666},
	}

	for _, c := range cases {
		got := c.rtsys.KillingForm(c.wt1, c.wt2)
		if got != c.want {
			t.Errorf("KillingForm(%v, %v) == %v, want %v", c.wt1, c.wt2, got, c.want)
		}
	}
}

func TestTypeAIntKillingForm(t *testing.T) {
	cases := []struct {
		rtsys    RootSystem
		wt1, wt2 Weight
		want     int
	}{
		{typeA{1}, Weight{0}, Weight{0}, 0},
		{typeA{1}, Weight{1}, Weight{0}, 0},
		{typeA{1}, Weight{0}, Weight{1}, 0},
		{typeA{1}, Weight{1}, Weight{1}, 1},
		{typeA{1}, Weight{2}, Weight{1}, 2},
		{typeA{1}, Weight{1}, Weight{2}, 2},
		{typeA{1}, Weight{2}, Weight{2}, 4},
		{typeA{2}, Weight{0, 0}, Weight{0, 0}, 0},
		{typeA{2}, Weight{1, 0}, Weight{0, 0}, 0},
		{typeA{2}, Weight{0, 0}, Weight{1, 0}, 0},
		{typeA{2}, Weight{1, 0}, Weight{1, 0}, 2},
		{typeA{2}, Weight{0, 1}, Weight{1, 0}, 1},
		{typeA{2}, Weight{1, 0}, Weight{0, 1}, 1},
		{typeA{2}, Weight{0, 1}, Weight{0, 1}, 2},
	}

	for _, c := range cases {
		got := c.rtsys.IntKillingForm(c.wt1, c.wt2)
		if got != c.want {
			t.Errorf("KillingForm(%v, %v) == %v, want %v", c.wt1, c.wt2, got, c.want)
		}
	}
}

func TestTypeAKillingFactor(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		want  int
	}{
		{typeA{1}, 2},
		{typeA{2}, 3},
		{typeA{3}, 4},
	}

	for _, c := range cases {
		got := c.rtsys.KillingFactor()
		if got != c.want {
			t.Errorf("DualCoxeter() == %v, want %v", got, c.want)
		}
	}
}

func TestTypeALevel(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		wt    Weight
		want  int
	}{
		{typeA{1}, Weight{0}, 0},
		{typeA{1}, Weight{1}, 1},
		{typeA{1}, Weight{2}, 2},
		{typeA{2}, Weight{0, 0}, 0},
		{typeA{2}, Weight{1, 0}, 1},
		{typeA{2}, Weight{0, 1}, 1},
		{typeA{2}, Weight{1, 1}, 2},
	}

	for _, c := range cases {
		got := c.rtsys.Level(c.wt)
		if got != c.want {
			t.Errorf("Level(%v) == %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeADual(t *testing.T) {
	cases := []struct {
		rtsys    RootSystem
		wt, want Weight
	}{
		{typeA{1}, Weight{0}, Weight{0}},
		{typeA{1}, Weight{1}, Weight{1}},
		{typeA{1}, Weight{2}, Weight{2}},
		{typeA{2}, Weight{0, 0}, Weight{0, 0}},
		{typeA{2}, Weight{1, 0}, Weight{0, 1}},
		{typeA{2}, Weight{0, 1}, Weight{1, 0}},
		{typeA{2}, Weight{1, 1}, Weight{1, 1}},
	}

	for _, c := range cases {
		got := c.rtsys.Dual(c.wt)
		if !equals(got, c.want) {
			t.Errorf("Dual(%v) = %v, want %v", c.wt, got, c.want)
		}
	}
}

func TestTypeAReflectIntoChamber(t *testing.T) {
	cases := []struct {
		rtsys    RootSystem
		wt, want Weight
		parity   int
	}{
		{typeA{1}, Weight{0}, Weight{0}, 1},
		{typeA{1}, Weight{1}, Weight{1}, 1},
		{typeA{1}, Weight{-1}, Weight{1}, -1},
		{typeA{2}, Weight{0, 0}, Weight{0, 0}, 1},
		{typeA{2}, Weight{1, 0}, Weight{1, 0}, 1},
		{typeA{2}, Weight{0, 1}, Weight{0, 1}, 1},
		{typeA{2}, Weight{-1, 0}, Weight{0, 1}, 1},
		{typeA{2}, Weight{0, -1}, Weight{1, 0}, 1},
		{typeA{2}, Weight{-1, -1}, Weight{1, 1}, -1},
	}

	for _, c := range cases {
		got := c.rtsys.NewWeight()
		parity := c.rtsys.reflectToChamber(c.wt, got)
		if !equals(got, c.want) || parity != c.parity {
			t.Errorf("ReflectToChamber(%v) = %v, %v, want %v, %v",
				c.wt, got, parity, c.want, c.parity)
		}
	}
}

func TestTypeAOrbitIterator(t *testing.T) {
	cases := []struct {
		rtsys RootSystem
		wt    Weight
		orbit []Weight
	}{
		{typeA{1}, Weight{0}, []Weight{Weight{0}}},
		{typeA{1}, Weight{1}, []Weight{Weight{1}, Weight{-1}}},
		{typeA{1}, Weight{2}, []Weight{Weight{2}, Weight{-2}}},
		{typeA{2}, Weight{0, 0}, []Weight{Weight{0, 0}}},
		{typeA{2}, Weight{1, 0}, []Weight{
			Weight{1, 0},
			Weight{-1, 1},
			Weight{0, -1}}},
		{typeA{2}, Weight{0, 1}, []Weight{
			Weight{0, 1},
			Weight{1, -1},
			Weight{-1, 0}}},
		{typeA{2}, Weight{1, 1}, []Weight{
			Weight{1, 1},
			Weight{-1, 2},
			Weight{2, -1},
			Weight{1, -2},
			Weight{-2, 1},
			Weight{-1, -1}}},
	}

	for _, c := range cases {
		orbitSet := weightSetFromList(c.orbit)
		var orbitEpc epCoord = make([]int, len(c.wt)+1)
		c.rtsys.convertWeightToEpc(c.wt, orbitEpc)
		orbitSize := 0
		done := false
		for ; !done; done = c.rtsys.nextOrbitEpc(orbitEpc) {
			nextWt := c.rtsys.NewWeight()
			c.rtsys.convertEpCoord(orbitEpc, nextWt)
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

func equals(v1, v2 []int) bool {
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
