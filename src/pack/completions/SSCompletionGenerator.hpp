/*
 * SchroeppelShamirCompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef SCHROEPPELSHAMIRCOMPLETIONGENERATOR_HPP_
#define SCHROEPPELSHAMIRCOMPLETIONGENERATOR_HPP_

#include "Completion_Generator.hpp"
#include "../../ss/SchroeppelShamir.hpp"
#include <deque>

using std::deque;
namespace ss {

// ===============================
// The Schroeppel and Shamir Class
// ===============================

class SSCompletionGenerator : public CompletionGenerator {
private :
  vector<SetNodeBitset> m_a;      // all subset sums from first quarter of numbers
  vector<SetNodeBitset> m_b;      // all subset sums from second quarter of numbers
  vector<SetNodeBitset> m_c;      // all subset sums from third quarter of numbers
  vector<SetNodeBitset> m_d;      // all subset sums from fourth quarter of numbers

  SSMinHeap m_abheap;        // min heap of sets from combination of A and B sets
  SSMaxHeap m_cdheap;        // max heap of sets from combination of C and D sets

  deque<SSHeapNode> m_cdlist; // list of subset sums removed from CD heap
  deque<SetNodeVector> m_generatedSets; // Sets that have been generated so far
  SSHeapNode m_newHeap;      // new element being added to heap

  uint64_t m_subsum;        // subset sum of first subset

  uint64_t m_lower;         // Lower Bound
  uint64_t m_upper;         // Upper Bound

  const int m_n;
  const uint64_t *m_S;

  bool m_generateAll;     // Generate all sets in one shot (true) or do it lazily (false)?


protected :
  bool generateNext();
public :

  SSCompletionGenerator(const uint64_t S[MAXN],const int n,
                        const uint64_t lower,const uint64_t upper,
                        bool generateAll = false);

  SSCompletionGenerator(const uint64_t S[MAXN],const int n,
                        bool generateAll = false);

  ~SSCompletionGenerator() {

  }

  bool next(SetNodeVector &node);

  void reset(const uint64_t lower,const uint64_t upper);
};

} /* namespace ss */
#endif /* SCHROEPPELSHAMIRCOMPLETION_HPP_ */
