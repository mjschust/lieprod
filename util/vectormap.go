package util

// VectorMap is a simple trie-based map from a float slice to arbitrary type values
type VectorMap interface {
	Get([]float64) (interface{}, bool)
	Put([]float64, interface{})
}

// NewVectorMap creates a new vector map.
func NewVectorMap() VectorMap {
	return &vectorMap{head: node{}}
}

type vectorMap struct {
	head node
}

func (vmap *vectorMap) Get(key []float64) (interface{}, bool) {
	curNode := &vmap.head
	for _, keyElem := range key {
		if curNode.children == nil {
			return nil, false
		} else {
			nextNode, present := curNode.children[keyElem]
			if present {
				curNode = nextNode
			} else {
				return nil, false
			}
		}
	}

	return curNode.value, curNode.value != nil
}

func (vmap *vectorMap) Put(key []float64, value interface{}) {
	curNode := &vmap.head
	for _, keyElem := range key {
		if curNode.children == nil {
			curNode.children = make(map[float64]*node)
			curNode.children[keyElem] = &node{}
			curNode = curNode.children[keyElem]
		} else {
			nextNode, present := curNode.children[keyElem]
			if present {
				curNode = nextNode
			} else {
				curNode.children[keyElem] = &node{}
				curNode = curNode.children[keyElem]
			}
		}
	}

	curNode.value = value
}

func (vmap *vectorMap) Remove(key []float64) {
	vmap.Put(key, nil)
}

type node struct {
	value    interface{}
	children map[float64]*node
}
