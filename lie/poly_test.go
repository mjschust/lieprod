package lie

import (
	"math/big"
	"testing"
)

func TestEmptyPoly(t *testing.T) {
	var emptyPoly WeightPoly = NewWeightPolyBuilder(1)
	if len(emptyPoly.Weights()) != 0 {
		t.Errorf("Empty weight poly contains weights: %v", emptyPoly.Weights())
	}

	wt := Weight{1}
	got := emptyPoly.Multiplicity(wt)
	if got.Cmp(big.NewInt(0)) != 0 {
		t.Errorf("Empty weight poly has non-zero multiplicity: emptyPoly[%v] = %v", wt, got)
	}
}

type monomial struct {
	wt           Weight
	multiplicity int
}

func (m monomial) mult() *big.Int {
	return big.NewInt(int64(m.multiplicity))
}

func TestSetMultiplicity(t *testing.T) {
	cases := []struct {
		rank int
		ms   []monomial
		want []monomial
	}{
		{
			1,
			[]monomial{},
			[]monomial{{Weight{1}, 0}},
		},
		{
			1,
			[]monomial{{Weight{1}, 1}},
			[]monomial{{Weight{1}, 1}},
		},
		{
			1,
			[]monomial{{Weight{1}, 1}, {Weight{2}, 2}},
			[]monomial{{Weight{1}, 1}, {Weight{2}, 2}},
		},
		{
			1,
			[]monomial{{Weight{1}, 1}, {Weight{1}, 2}},
			[]monomial{{Weight{1}, 2}},
		},
	}

	for _, c := range cases {
		poly := NewWeightPolyBuilder(c.rank)
		for _, mono := range c.ms {
			poly.SetMultiplicity(mono.wt, mono.mult())
		}

		for _, mono := range c.want {
			gotMult := poly.Multiplicity(mono.wt)
			if gotMult.Cmp(mono.mult()) != 0 {
				t.Errorf("poly[%v] = %v, want %v", mono.wt, gotMult, mono.mult())
			}
		}
	}
}

func TestAddMultiplicity(t *testing.T) {
	cases := []struct {
		rank int
		ms   []monomial
		want []monomial
	}{
		{
			1,
			[]monomial{{Weight{1}, 1}},
			[]monomial{{Weight{1}, 1}},
		},
		{
			1,
			[]monomial{{Weight{1}, 1}, {Weight{2}, 2}},
			[]monomial{{Weight{1}, 1}, {Weight{2}, 2}},
		},
		{
			1,
			[]monomial{{Weight{1}, 1}, {Weight{1}, 2}},
			[]monomial{{Weight{1}, 3}},
		},
		{
			2,
			[]monomial{{Weight{1, 1}, 1}, {Weight{1, 1}, 2}},
			[]monomial{{Weight{1, 1}, 3}},
		},
	}

	for _, c := range cases {
		poly := NewWeightPolyBuilder(c.rank)
		for _, mono := range c.ms {
			poly.AddMultiplicity(mono.wt, mono.mult())
		}

		for _, mono := range c.want {
			gotMult := poly.Multiplicity(mono.wt)
			if gotMult.Cmp(mono.mult()) != 0 {
				t.Errorf("poly[%v] = %v, want %v", mono.wt, gotMult, mono.mult())
			}
		}
	}
}

func TestSingleWeightPoly(t *testing.T) {
	poly := Weight{1}

	gotWts := poly.Weights()
	if len(gotWts) != 1 || !gotWts[0].Equals(Weight{1}) {
		t.Errorf("For %v, WeightPoly.GetWeights() = %v, want %v", poly, gotWts, []Weight{{1}})
	}

	gotMult := poly.Multiplicity(Weight{1})
	if gotMult.Cmp(big.NewInt(1)) != 0 {
		t.Errorf("For %v, WeightPoly.GetMultiplicity(%v) = %v, want %v",
			poly, Weight{1}, gotMult, 1)
	}

	gotMult = poly.Multiplicity(Weight{0})
	if gotMult.Cmp(big.NewInt(0)) != 0 {
		t.Errorf("For %v, WeightPoly.GetMultiplicity(%v) = %v, want %v",
			poly, Weight{0}, gotMult, 0)
	}

	gotMult = poly.Multiplicity(Weight{1, 2})
	if gotMult.Cmp(big.NewInt(0)) != 0 {
		t.Errorf("For %v, WeightPoly.GetMultiplicity(%v) = %v, want %v",
			poly, Weight{1, 2}, gotMult, 0)
	}
}
