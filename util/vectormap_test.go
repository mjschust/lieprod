package util

import (
	"testing"

	"github.com/mjschust/cblocks/lie"
)

func TestEmptyEntries(t *testing.T) {
	cases := []struct {
		existingKey []float64
		testKey     []float64
	}{
		{[]float64{}, []float64{1}},
		{[]float64{1}, []float64{}},
		{[]float64{1}, []float64{1, 2}},
		{[]float64{1, 2}, []float64{}},
		{[]float64{1, 2}, []float64{1}},
		{[]float64{1, 2}, []float64{1, 2, 3}},
	}

	for _, c := range cases {
		vmap := NewVectorMap()
		vmap.Put(c.existingKey, 1)
		got, _ := vmap.Get(c.testKey)
		if got != nil {
			t.Errorf("vmap.Get(%v) = %v, want nil", c.testKey, got)
		}
	}
}
func TestSinglePut(t *testing.T) {
	cases := []struct {
		key   []float64
		value interface{}
	}{
		{[]float64{1}, 5},
		{lie.Weight{2}, true},
		{[]float64{1, 0}, 5},
		{[]float64{1, 1}, 1.5},
	}

	for _, c := range cases {
		vmap := NewVectorMap()
		vmap.Put(c.key, c.value)
		got, _ := vmap.Get(c.key)
		if got != c.value {
			t.Errorf("vmap.Get(%v) = %v, want %v", c.key, got, c.value)
		}
	}
}

func TestMultiPut(t *testing.T) {
	cases := []struct {
		key   []float64
		value interface{}
	}{
		{[]float64{1}, 1},
		{[]float64{1, 0}, 2},
		{[]float64{1, 1}, 2.5},
		{[]float64{2, 1}, 3.5},
	}

	vmap := NewVectorMap()
	for _, c := range cases {
		vmap.Put(c.key, c.value)
		got, _ := vmap.Get(c.key)
		if got != c.value {
			t.Errorf("vmap.Get(%v) = %v, want %v", c.key, got, c.value)
		}
	}
}
