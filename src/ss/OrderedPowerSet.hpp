/*
 * OrderedPowerSet.hpp
 *
 *  Created on: Jun 16, 2013
 *      Author: ethan
 */

#ifndef ORDEREDPOWERSET_HPP_
#define ORDEREDPOWERSET_HPP_
#include "../Heap.hpp"
#include <stdint.h>
#include <sstream>
#include <vector>
#include <limits>

namespace ss {

// =======================================================
// The Node to be stored in the OrderedPowerSet Heap class
// =======================================================
struct OPSHeapNode {
  size_t LIdx;        // index of set in L
  size_t SIdx;        // index of element in S
  uint64_t sum;      // sum of the elements in subset + element
  int cardinality;   // cardinality of set of elements

  OPSHeapNode(size_t LIdx0, size_t SIdx0, uint64_t sum0, int cardinality0) :
    LIdx(LIdx0), SIdx(SIdx0), sum(sum0), cardinality(cardinality0) {}

  OPSHeapNode(const OPSHeapNode &node) :
    LIdx(node.LIdx), SIdx(node.SIdx), sum(node.sum), cardinality(node.cardinality) {}
};

// ====================================
// Min and Max heap Comparators for OPS
// ====================================
struct OPSMinHeapComparator {

public :
  bool operator()(const OPSHeapNode &a, const OPSHeapNode &b) {
    if (a.sum == b.sum) {   // if sums equal, return smaller cardinality first
      return a.cardinality > b.cardinality;
    } else {
      return a.sum > b.sum;
    }
  }
};

struct OPSMaxHeapComparator {

public :
  bool operator()(const OPSHeapNode &a, const OPSHeapNode &b) {
    if (a.sum == b.sum) { // if sums equal, return smaller cardinality first
      return a.cardinality > b.cardinality;
    } else {
      return a.sum < b.sum;
    }
  }
};

// ============================================
// Typedef MinHeap and MaxHeap to use the heaps
// ============================================
typedef Heap<OPSMinHeapComparator,OPSHeapNode> OPSMinHeap;
typedef Heap<OPSMaxHeapComparator,OPSHeapNode> OPSMaxHeap;

//
template <typename Type>
class PersisistentIndexDeque {
private :
  deque<Type> &m_d;
  size_t m_offset;
  size_t m_maxCount;

public :
  PersisistentIndexDeque(deque<SetNodeBitset> &d) : m_d(d), m_offset(0), m_maxCount(0) {
    m_d.clear();
  }

  void push_back(Type value) {
    m_d.push_back(value);
    m_maxCount = std::max(m_maxCount,m_d.size());
  }

  size_t maxCount() const {
  	return m_maxCount;
  }

  Type& back() const {
    return m_d.back();
  }

  void pop_front() {
    m_offset++;
    m_d.pop_front();

  }

  inline Type& operator[](size_t idx){
    return m_d[idx - m_offset];
  }


  size_t pushCount() const {
    return m_d.size() + m_offset;
  }

  size_t count() const {
    return m_d.size();
  }
};

}

#endif /* ORDEREDPOWERSET_HPP_ */
