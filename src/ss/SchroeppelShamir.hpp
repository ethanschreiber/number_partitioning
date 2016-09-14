/*
 * Schroeppel_Shamir.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: ethan
 */

#ifndef SCHROEPPEL_SHAMIR_HPP_
#define SCHROEPPEL_SHAMIR_HPP_

#include "SubsetSum.hpp"
#include "../partition/PartitionUtils.hpp"
#include "../Heap.hpp"
#include <vector>
#include <deque>
#include <algorithm>

using std::vector;
using std::deque;

namespace ss {

// ============================================================
// The Node to be stored in the Schroeppel and Shmir Heap class
// ============================================================
class SSHeapNode {   // heap of sums from a and b in increasing order
private :
  size_t m_x;           // index of element in asums
  size_t m_y;           // index of element in bsums
  uint64_t m_sum;        // sum of the elements

public :
  SSHeapNode(size_t x0, size_t y0,uint64_t sum0) : m_x(x0), m_y(y0), m_sum(sum0) {}

  SSHeapNode(size_t x0, size_t y0,const ss::SetNodeBitset &setNode) : m_x(x0), m_y(y0), m_sum(setNode.sum()) {}

  SSHeapNode(size_t x0, size_t y0,const vector<uint64_t> &xVector,const vector<uint64_t> &yVector)
  : m_x(x0), m_y(y0), m_sum(xVector[x0] + yVector[y0]) { }

  SSHeapNode(size_t x0, size_t y0,const vector<ss::SetNodeBitset> &xVector,const vector<ss::SetNodeBitset> &yVector)
  : m_x(x0), m_y(y0), m_sum(xVector[x0].sum() + yVector[y0].sum()) { }

  SSHeapNode() {}

  void set(size_t x, size_t y, uint64_t sum) {
    m_x = x;
    m_y = y;
    m_sum = sum;
  }
  const size_t& x() const{
    return m_x;
  }

  const size_t& y() const{
    return m_y;
  }

  const size_t& sum() const{
    return m_sum;
  }
};

// ======================================================
// Min and Max heap Comparators for Schroeppel and Shamir
// ======================================================

template <typename NodeTypeT>
class SSMinHeapComparator {

public :
  bool operator()(const NodeTypeT &a, const NodeTypeT &b) {
    return a.sum() > b.sum();
  }
};

template <typename NodeTypeT>
class SSMaxHeapComparator {

public :
  bool operator()(const NodeTypeT &a, const NodeTypeT &b) {
    return a.sum() < b.sum();
  }
};

// ============================================
// Typedef MinHeap and MaxHeap to use the heaps
// ============================================
typedef Heap<SSMinHeapComparator<SSHeapNode>,SSHeapNode> SSMinHeap;
typedef Heap<SSMaxHeapComparator<SSHeapNode>,SSHeapNode> SSMaxHeap;

typedef Heap<SSMinHeapComparator<SetNodeBitset>,SetNodeBitset> SSSetMinHeap;
typedef Heap<SSMaxHeapComparator<SetNodeBitset>,SetNodeBitset> SSSetMaxHeap;


// ============================================================================
// Schroeppel and Shamir Helper Functions for *** SetNodeBitset ***
// ============================================================================

void initializeSets(const uint64_t S[MAXN], const int n,
                    vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                    const uint8_t maxCardinality= UNSET_MAX_CARD);

void initializeHeaps(vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                     SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper, const uint8_t maxCardinality= UNSET_MAX_CARD);

void initialize(const uint64_t S[MAXN], const int n,
                       vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                       SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper, const uint8_t maxCardinality= UNSET_MAX_CARD);

// ============================================================================
// Schroeppel and Shamir Helper Functions for *** Sums Only ***
// ============================================================================

void initializeSets(const uint64_t S[MAXN], const int n,
                           vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d);

void initializeHeaps(vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d,
                     SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper);

void initialize(const uint64_t S[MAXN], const int n,
                       vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d,
                       SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper);


// ============================================================================
// Schroeppel and Shamir Partition Problem
// ============================================================================

/**
 * Run Schroeppel and Shamir to solve the (two-way) Partition problem.
 */
uint64_t executeSS2(const partition::PartitionProblem &problem, ProblemStats &stats);


} // end namespace



#endif /* SCHROEPPEL_SHAMIR_HPP_ */
