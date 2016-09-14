/*
 * VectorCompletion.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#include "OPSCompletionGenerator.hpp"


namespace ss {

OPSCompletionGenerator::OPSCompletionGenerator(const uint64_t *S, const int N,
    const uint64_t lower, const uint64_t upper, deque<SetNodeBitset> &LParam)
: m_N(N),m_S(S), m_lower(lower), m_upper(upper), m_L(LParam) {

  m_L.push_back(SetNodeBitset(N));                   // Add empty subset

  // Initialize the heap
  for (int i=0;i<N;i++) {                     // For each input element
    m_nextHeap.push(OPSHeapNode(0,i,S[i],0)); // One entry in heap pointing to 0 idx of m_L, i idx of S, sum of S[i] and cardinality 0
  }
}

OPSCompletionGenerator::~OPSCompletionGenerator() {
}

bool OPSCompletionGenerator::next(SetNodeVector &node) {
  OPSHeapNode heapNode = m_nextHeap.front();            // Get min sum
  m_nextHeap.pop();                                     // Pop it from heap

  if (heapNode.sum <= m_upper &&                        // Add if less than m_upper bound AND
      (heapNode.sum >= m_lower ||                       // greater than m_upper bound OR
       heapNode.SIdx < (size_t) m_N-1)) {               // another set can be created from this one.
    m_L.push_back(SetNodeBitset(m_L[heapNode.LIdx]));          // 3. Insert copy of old entry into m_L
    m_L.back().set(heapNode.SIdx, m_S[heapNode.SIdx]);  // Insert new element into copy (changes sum too)
  }

  if (heapNode.SIdx == (size_t) m_N-1 &&                // If this node corresponds to the largest value in S
      m_L[heapNode.LIdx].sum() < m_lower) {             // And the value is less than the m_lower bound

    m_L.pop_front();
  }

  size_t LIdx = heapNode.LIdx;
  do {                                                  // Find next set we can add S[node.SIdx] to.
    LIdx++;
  } while (LIdx < m_L.pushCount() && m_L[LIdx].max() >= heapNode.SIdx);

  if (LIdx < m_L.pushCount()) {      // If there is another subset to point to
    m_nextHeap.push(OPSHeapNode(LIdx,heapNode.SIdx,m_L[LIdx].sum() + m_S[heapNode.SIdx],m_L[LIdx].set().count()));
  }

  return !m_nextHeap.empty();
}



} /* namespace ss */
