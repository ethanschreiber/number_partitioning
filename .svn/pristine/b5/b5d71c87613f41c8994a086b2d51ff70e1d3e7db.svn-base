/*
 * SchroeppelShamirCompletion.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#include "SSCompletionGenerator.hpp"
#include "../BinCompletionUtils.hpp"
#include "../PackingUtils.hpp"
#include <algorithm>
namespace ss {

// ===============================
// The Schroeppel and Shamir Class
// ===============================

SSCompletionGenerator::SSCompletionGenerator(
    const uint64_t S[MAXN],const int n,const uint64_t lower,const uint64_t upper, bool generateAll)
: m_lower(lower), m_upper(upper), m_n(n), m_S(S),m_generateAll(generateAll)
{
  initialize(S,n,m_a,m_b,m_c,m_d,m_abheap,m_cdheap,upper);
}

SSCompletionGenerator::SSCompletionGenerator(
    const uint64_t S[MAXN],const int n, bool generateAll) : m_n(n), m_S(S),m_generateAll(generateAll)
{
  initializeSets(m_S,m_n,m_a,m_b,m_c,m_d);
}

void SSCompletionGenerator::reset(const uint64_t lower,const uint64_t upper) {
  m_lower = lower;
  m_upper = upper;
  m_abheap.clear();
  m_cdheap.clear();
  m_cdlist.clear();
  initializeHeaps(m_a,m_b,m_c,m_d,m_abheap,m_cdheap,m_upper);
}

bool SSCompletionGenerator::generateNext() {
  bool foundAny = false;
  while (!m_abheap.empty()) {                      // until AB heap is empty

    // remove sets with sums above upper bound from m_cdlist
    while (!m_cdlist.empty() && m_abheap.front().sum() + m_cdlist.front().sum() > m_upper) {
      m_cdlist.pop_front();
    }

    // add new sets with sums within bounds to m_cdlist
    uint64_t topsum = m_cdheap.front().sum() + m_abheap.front().sum();     //sum of elements on top of heaps
    while (!m_cdheap.empty() && topsum >= m_lower) {       	// new set within bounds

      if (topsum <= m_upper) {                             	// new set within bounds
        m_cdlist.push_back(m_cdheap.front());              	// add heap element to list
      }

      if (m_cdheap.front().y() > 0) {                        	// there is another element in this column
        SSHeapNode newElement(m_cdheap.front().x(),          	// keep same row
                            m_cdheap.front().y() - 1,        	// replace with next combination in column
                            m_c,m_d);                      	// sum is c[x] + d[y]

        m_cdheap.pop();                                    	// Remove old head
        m_cdheap.push(newElement);                         	// add new element to heap vector

      } else {                                           	// no more elements in this column
        m_cdheap.pop();                                    	// Remove old head
      }

      topsum = m_cdheap.front().sum() + m_abheap.front().sum();      // sum of two top elements of heaps
    }

    // Combine each element from m_cdlist with the top element from
    // m_abheap, these are all within range
    for (size_t index=0;index<m_cdlist.size();index++) {
      SetNodeVector included;
      SetNodeVector excluded;

      DynamicBitset bitset = m_c[m_cdlist[index].x()].set() |
                                    m_d[m_cdlist[index].y()].set() |
                                    m_a[m_abheap.front().x()].set() |
                                    m_b[m_abheap.front().y()].set();

      int leastSetIdx = m_n-1;            // Find the least significant bit set

      while (!bitset[leastSetIdx]) {  // Assume at least 1 bit is set
        leastSetIdx--;
      }

      for (int idx = 0;idx<=leastSetIdx;idx++) {
        if (bitset[idx]) {
          included.push_back(m_S[idx]);
        } else {
          excluded.push_back(m_S[idx]);
        }
      }

      if (!bp::isDominated(included.elements(), included.getSum(), excluded.elements(),m_upper-included.getSum())) {

        m_generatedSets.push_back(included);
        foundAny = true;
      }
    }

    // If there is another element from the b set, add it.
    // Then, remove top element from m_abheap.
    if (m_abheap.front().y() < m_b.size() - 1) {          // there is another element in this column
      SSHeapNode newNode(m_abheap.front().x(),            // new element being added to heap
                       m_abheap.front().y() + 1,        // replace with next combination in this column
                       m_a,m_b);

      m_abheap.pop();                            // Pop old value
      m_abheap.push(newNode);                    // add new element to heap vector
    } else {
      m_abheap.pop();
    }
    if (foundAny && !m_generateAll) {
      break;
    }
  }
  return foundAny;
}


bool SSCompletionGenerator::next(SetNodeVector &node) {

  bool foundAny = !m_generatedSets.empty();
  if (!foundAny) {
    foundAny = generateNext();
    std::sort(m_generatedSets.begin(),m_generatedSets.end(),SetNodeCardinalityComparator());
  }

  if (foundAny) {                             // If any sets
//    node.set = m_generatedSets.front().set;   // Copy them
//    node.sum = m_generatedSets.front().sum;

    node = m_generatedSets.front();
    m_generatedSets.pop_front();              // Remove it
  }

  return foundAny;
}


} /* namespace ss */
