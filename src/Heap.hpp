/*
 * Heap.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef HEAP_HPP
#define HEAP_HPP

#include "ss/SubsetSum.hpp"
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

using std::vector;
using std::binary_function;



// ===================
// Generic number heap
// ===================

template <class T>
class MinHeapComparator {
public:
	bool operator() ( const T &a, const T &b ) const {
		return a > b;
	}
};

template <class T>
class MaxHeapComparator {
public:
	bool operator() ( const T &a, const T &b ) const {
		return a < b;
	}
};


// =====================================================
// The Heap data structure has a vector underlying
// Add data initially using the vector interface.
// Then, when ready, call make_heap() and then only use
// push_heap and pop_heap
// =====================================================
// Class C is the comparator
template <class Comparator, typename Node>
class Heap {
private :
  Comparator  m_comparator;
  vector<Node> m_vector;

public:

  void make_heap() {
    std::make_heap(m_vector.begin(),m_vector.end(),m_comparator);
  }

  void pop() {
    std::pop_heap (m_vector.begin(),m_vector.end(),m_comparator); // Remove old head
    m_vector.pop_back();                                 // remove it from vector
  }

  void push(const Node &node) {
  	m_vector.push_back(node);                             // add new element to heap vector
    std::push_heap(m_vector.begin(),m_vector.end(),m_comparator);  // percolate to proper position
  }

  Node &front() {
	  return m_vector.front();
  }

  // This just adds elements to the vector, use this before calling make_heap
  void vector_push_back (const Node& node) {
	  m_vector.push_back(node);
  }
  bool empty() const {
	  return m_vector.empty();
  }

  void clear() {
  	m_vector.clear();
  }

  size_t size() const{
	  return m_vector.size();
  }

  virtual ~Heap() {}
};

// ============================================
// Typedef MinHeap and MaxHeap to use the heaps
// ============================================

typedef Heap<MinHeapComparator<int>,int> IntMinHeap;
typedef Heap<MaxHeapComparator<int>,int> IntMaxHeap;

#endif /* HEAP_HPP_ */
