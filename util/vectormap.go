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
	return &vectorMap{head: &node{}}
}

type vectorMap struct {
	head *node
	size int
}

func (vmap *vectorMap) Get(key []int) (interface{}, bool) {
	curNode := vmap.head
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

func (vmap *vectorMap) Put(key []int, value interface{}) {
	if value == nil {
		return
	}

	newNode := false
	curNode := vmap.head
	for _, keyElem := range key {
		if curNode.children == nil {
			curNode.children = make(map[int]*node)
			curNode.children[keyElem] = &node{key: keyElem, parent: curNode}
			curNode = curNode.children[keyElem]
			newNode = true
		} else {
			nextNode, present := curNode.children[keyElem]
			if present {
				curNode = nextNode
			} else {
				curNode.children[keyElem] = &node{key: keyElem, parent: curNode}
				curNode = curNode.children[keyElem]
				newNode = true
			}
		}
	}

	if newNode || curNode.value == nil {
		vmap.size++
	}
	curNode.value = value
}

func (vmap *vectorMap) Remove(key []int) (interface{}, bool) {
	curNode := vmap.head
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

	oldValue := curNode.value
	curNode.value = nil
	vmap.removeEmptyNode(curNode)

	if oldValue != nil {
		vmap.size--
	}

	return oldValue, oldValue != nil
}

func (vmap *vectorMap) removeEmptyNode(n *node) {
	switch {
	case n.value != nil:
		return
	case n.parent == nil:
		return
	case n.children != nil && len(n.children) > 0:
		return
	default:
		delete(n.parent.children, n.key)
		vmap.removeEmptyNode(n.parent)
	}
}

func (vmap *vectorMap) Keys() [][]int {
	keys := [][]int{}

	var keysHelper func(n *node, keyPrefix []int)
	keysHelper = func(n *node, keyPrefix []int) {
		if n.children != nil && len(n.children) > 0 {
			newPrefix := make([]int, len(keyPrefix)+1)
			copy(newPrefix, keyPrefix)
			for keyElem, childNode := range n.children {
				newPrefix[len(keyPrefix)] = keyElem
				if childNode.value != nil {
					key := make([]int, len(newPrefix))
					copy(key, newPrefix)
					keys = append(keys, key)
				}
				keysHelper(childNode, newPrefix)
			}
		}
	}

	keysHelper(vmap.head, []int{})
	return keys
}

func (vmap *vectorMap) Size() int {
	return vmap.size
}

type node struct {
	key      int
	value    interface{}
	parent   *node
	children map[int]*node
}
