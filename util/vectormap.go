package util

// VectorMap is a simple trie-based map from a float slice to arbitrary type values
type VectorMap interface {
	Get([]int) (interface{}, bool)
	Put([]int, interface{})
	Remove([]int) (interface{}, bool)
	Keys() [][]int
	Size() int
}

// NewVectorMap creates a new vector map.
func NewVectorMap() VectorMap {
	return &vectorMap{table: make(map[int][]entry)}
}

type vectorMap struct {
	table map[int][]entry
	size  int
}

func (vmap *vectorMap) Get(key []int) (interface{}, bool) {
	hash := hashKey(key)
	bucket := vmap.table[hash]
	if bucket != nil {
		_, e := findEntry(key, bucket)
		if e != nil {
			return e.getValue(), e.getValue() != nil
		}
	}

	return nil, false
}

func (vmap *vectorMap) Put(key []int, value interface{}) {
	hash := hashKey(key)
	bucket := vmap.table[hash]
	if bucket != nil {
		_, e := findEntry(key, bucket)
		if e != nil {
			e.setValue(value)
		} else {
			vmap.table[hash] = append(bucket, &mapEntry{key, value})
			vmap.size++
		}
	} else {
		vmap.table[hash] = []entry{&mapEntry{key, value}}
		vmap.size++
	}
}

func (vmap *vectorMap) Remove(key []int) (interface{}, bool) {
	hash := hashKey(key)
	bucket := vmap.table[hash]
	if bucket != nil {
		i, e := findEntry(key, bucket)
		if e != nil {
			vmap.table[hash] = removeEntry(bucket, i)
			vmap.size--
			return e.getValue(), e.getValue() != nil
		}
	}

	return nil, false
}

func (vmap *vectorMap) Keys() [][]int {
	keys := make([][]int, 0, vmap.size)
	for _, bucket := range vmap.table {
		for _, e := range bucket {
			keys = append(keys, e.getKey())
		}
	}
	return keys
}

func (vmap *vectorMap) Size() int {
	return vmap.size
}

func findEntry(key []int, bucket []entry) (int, entry) {
	for i, ent := range bucket {
		if equals(key, ent.getKey()) {
			return i, ent
		}
	}
	return 0, nil
}

func hashKey(key []int) int {
	retVal := 1
	for _, elem := range key {
		retVal = retVal*31 + elem
	}
	return retVal
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

func removeEntry(slice []entry, s int) []entry {
	return append(slice[:s], slice[s+1:]...)
}

type entry interface {
	getKey() []int
	getValue() interface{}
	setValue(interface{})
}

type mapEntry struct {
	key   []int
	value interface{}
}

func (e *mapEntry) getKey() []int {
	return e.key
}

func (e *mapEntry) getValue() interface{} {
	return e.value
}

func (e *mapEntry) setValue(value interface{}) {
	e.value = value
}
