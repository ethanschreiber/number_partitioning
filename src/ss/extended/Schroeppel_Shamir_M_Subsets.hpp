/*
 * Schroeppel_Shamir.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: ethan
 */

#ifndef SCHROEPPEL_SHAMIR_M_SUBSETS_HPP_
#define SCHROEPPEL_SHAMIR_M_SUBSETS_HPP_

#include "../SchroeppelShamir.hpp"

namespace ss {

// This one for Cached Iterative Weakening
// Generate the $m$ smallest greater than perfect

// ============================================
// Typedef MinHeap and MaxHeap to use the heaps
// ============================================


  class SchroeppelShamirMSubsets {
  private :

    const int m_K;
    // We split our N input elements into four quartiles. We generate
    // all 2^(n/4) subsets from each set and put them into
    // m_a, m_b, m_c and m_d each in sorted order
    vector<SetNodeBitset> m_a;  // The s^{n/4} subsets from the smallest quartile
    vector<SetNodeBitset> m_b;  // The s^{n/4} subsets from the second quartile
    vector<SetNodeBitset> m_c;  // The s^{n/4} subsets from the third quartile
    vector<SetNodeBitset> m_d;  // The s^{n/4} subsets from the largest quartile

    size_t m_numSets;     			// The number of sets to generate
    size_t m_elementsSum;				// The sum of the elements of m_S
    size_t m_lowerBound;  			// The current lower bound and
    size_t m_upperBound;  			// upper bound for generating sums

    size_t m_perfectValue; 			// The value of a perfect partition
    uint64_t m_smallMax;   			// The max subset sum of the subsets < perfect
    uint64_t m_largeMin;   			// The min subset sum of the subsets > perfect

    bool m_firstCall;     			// Double m_numSets after the first call
    double m_memoryUsage;				// Keep track of memory usage.
    double m_maxCard;           // The max cardinality allowed for a set
  public :

    // If you don't pass maxCard, then the max card is unlimited
    // (In reality, the max card is still sizeof(uint8_t = 255)
    // but Schroeppel and Shamir doesn't work memory wise with
    // more than about 100 ints)
    SchroeppelShamirMSubsets(const uint64_t S[], const int N, const int K, const size_t numSets, uint64_t upperBound, uint8_t maxCard=UNSET_MAX_CARD);

    // smallSets and largeSets are filled in.
    // smallSets contain all of the sets generated that are smaller than perfect
    //    - ordered from largest sum to smallest sum
    // largeSets contain all of the sets generated that are larger than perfect
    //    - ordered from smallest sum to largest sum
    void  generateSetsSS(vector<ss::SetNodeBitset> &smallSets, vector<ss::SetNodeBitset> &largeSets);

    inline double getMemoryUsage() const {
    	return m_memoryUsage;
    }
  };

} // end namespace

#endif /* SCHROEPPEL_SHAMIR_CARDINALITY_HPP_ */
