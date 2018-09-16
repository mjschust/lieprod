package util

import (
	"testing"
)

func TestEmptyEntries(t *testing.T) {
	cases := []struct {
		existingKey []int
		testKey     []int
	}{
		{[]int{}, []int{1}},
		{[]int{1}, []int{}},
		{[]int{1}, []int{1, 2}},
		{[]int{1, 2}, []int{}},
		{[]int{1, 2}, []int{1}},
		{[]int{1, 2}, []int{1, 2, 3}},
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
		key   []int
		value interface{}
	}{
		{[]int{1}, 5},
		{[]int{2}, true},
		{[]int{1, 0}, 5},
		{[]int{1, 1}, 1.5},
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
		key   []int
		value interface{}
	}{
		{[]int{1}, 1},
		{[]int{1, 0}, 2},
		{[]int{1, 1}, 2.5},
		{[]int{2, 1}, 3.5},
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

func TestRemove(t *testing.T) {
	cases := []struct {
		keys         [][]int
		removeKeys   [][]int
		wantPresent  [][]int
		wantAbsent   [][]int
		wantNumNodes int
	}{
		{[][]int{{1}}, [][]int{{1}}, [][]int{}, [][]int{{1}}, 0},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2, 3}}, [][]int{{1}, {1, 2}}, [][]int{{1, 2, 3}}, 2},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}, {1, 2, 3}}, [][]int{{1}}, [][]int{{1, 2}, {1, 2, 3}}, 1},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{}, [][]int{{1}, {1, 2}, {1, 2, 3}}, 0},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}}, [][]int{{1}, {1, 2, 3}}, [][]int{{1, 2}}, 3},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}}, [][]int{{1, 2}, {1, 2, 3}}, [][]int{{1}}, 3},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 3}}, [][]int{{1}, {1, 2}}, [][]int{{1, 3}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}}, [][]int{{1}, {1, 3}}, [][]int{{1, 2}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}, {1, 3}}, [][]int{{1}}, [][]int{{1, 2}, {1, 3}}, 1},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}}, [][]int{{1, 3}, {1, 2}}, [][]int{{1}}, 3},
	}

	for _, c := range cases {
		vmap := NewVectorMap()
		for i := range c.keys {
			vmap.Put(c.keys[i], true)
		}

		for _, key := range c.removeKeys {
			_, present := vmap.Remove(key)
			if !present {
				t.Errorf("vmap.Remove(%v) = _, false, want _, true", key)
			}
		}

		for _, key := range c.wantPresent {
			_, present := vmap.Get(key)
			if !present {
				t.Errorf("After removal, vmap.Get(%v) = _, false, want _, true", key)
			}
		}

		for _, key := range c.wantAbsent {
			_, present := vmap.Get(key)
			if present {
				t.Errorf("After removal, vmap.Get(%v) = _, true, want _, false", key)
			}
		}
	}
}

func TestKeys(t *testing.T) {
	cases := []struct {
		keys       [][]int
		removeKeys [][]int
		wantKeys   [][]int
	}{
		{[][]int{{1}}, [][]int{}, [][]int{{1}}},
		{[][]int{{1}}, [][]int{{1}}, [][]int{}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{}, [][]int{{1}, {1, 2}, {1, 2, 3}}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2, 3}}, [][]int{{1}, {1, 2}}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}, {1, 2, 3}}, [][]int{{1}}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}}, [][]int{{1}, {1, 2, 3}}},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}}, [][]int{{1, 2}, {1, 2, 3}}},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{}, [][]int{{1}, {1, 2}, {1, 3}}},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 3}}, [][]int{{1}, {1, 2}}},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}}, [][]int{{1}, {1, 3}}},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}, {1, 3}}, [][]int{{1}}},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}}, [][]int{{1, 3}, {1, 2}}},
	}

	for _, c := range cases {
		vmap := NewVectorMap()
		for i := range c.keys {
			vmap.Put(c.keys[i], true)
		}

		for _, key := range c.removeKeys {
			_, present := vmap.Remove(key)
			if !present {
				t.Errorf("vmap.Remove(%v) = _, false, want _, true", key)
			}
		}

		// Test that size of Keys() is equal to Size()
		if len(vmap.Keys()) != vmap.Size() {
			t.Errorf("%v : vmap.Keys() is size %v, want %v", c, len(vmap.Keys()), vmap.Size())
		}

		// Test that all keys in Keys() are in the map
		newVmap := NewVectorMap()
		for _, key := range vmap.Keys() {
			_, present := vmap.Get(key)
			if !present {
				t.Errorf("For %v in Keys(): vmap.Get(%v) = _, false, want _, true", key, key)
			}
			newVmap.Put(key, true)
		}

		// Test that all expected keys are in Keys()
		for _, key := range c.wantKeys {
			_, present := newVmap.Get(key)
			if !present {
				t.Errorf("vmap.Keys() is missing key %v", key)
			}
		}
	}
}

func TestSize(t *testing.T) {
	cases := []struct {
		keys       [][]int
		removeKeys [][]int
		want       int
	}{
		{[][]int{}, [][]int{}, 0},
		{[][]int{{1}}, [][]int{}, 1},
		{[][]int{{1}}, [][]int{{1}}, 0},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{}, 3},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2, 3}}, 2},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}, {1, 2, 3}}, 1},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}, {1, 2}, {1, 2, 3}}, 0},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1, 2}}, 2},
		{[][]int{{1}, {1, 2}, {1, 2, 3}}, [][]int{{1}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{}, 3},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 3}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1, 2}, {1, 3}}, 1},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}}, 2},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}, {1, 3}}, 1},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}, {1, 2}}, 1},
		{[][]int{{1}, {1, 2}, {1, 3}}, [][]int{{1}, {1, 2}, {1, 3}}, 0},
	}

	for _, c := range cases {
		vmap := NewVectorMap()
		for i := range c.keys {
			vmap.Put(c.keys[i], true)
		}

		for _, key := range c.removeKeys {
			_, present := vmap.Remove(key)
			if !present {
				t.Errorf("vmap.Remove(%v) = _, false, want _, true", key)
			}
		}

		got := vmap.Size()
		if got != c.want {
			t.Errorf("%v : Size() = %v, want %v", c, got, c.want)
		}
	}
}

func BenchmarkPut(b *testing.B) {
	vecs := make([][]int, 1)
	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			for k := 0; k < 10; k++ {
				for l := 0; l < 10; l++ {
					vecs = append(vecs, []int{i, j, k, l})
				}
			}
		}
	}

	b.ReportAllocs()
	b.ResetTimer()
	vmap := NewVectorMap()
	for i := 0; i < b.N; i++ {
		for i, vec := range vecs {
			vmap.Put(vec, i)
		}
	}
}
