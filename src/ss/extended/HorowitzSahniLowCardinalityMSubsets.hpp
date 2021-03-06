/*
 * Schroeppel_Shamir.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: ethan
 */

#ifndef HOROWITZ_SAHNI_LOW_CARDINALIITY_M_SUBSETS
#define HOROWITZ_SAHNI_LOW_CARDINALIITY_M_SUBSETS

#include "../SubsetSum.hpp"

namespace ss {

// This one for Low Cardinality Cached Iterative Weakening
// Generate the $m$ smallest greater than perfect with cardinality less
// than a MAX_CARDINALITY

// ============================================
// Typedef MinHeap and MaxHeap to use the heaps
// ============================================

  class HorowitzSahniLowCardinalityMSubsets {
  private :
    const std::vector<uint64_t> m_S;
    const int m_K;
    // We split our N input elements into four quartiles. We generate
    // all 2^(n/4) subsets from each set and put them into
    // m_a, m_b, m_c and m_d each in sorted order
    vector<SetNodeBitset> m_a;  // The subsets from the lower half
    vector<SetNodeBitset> m_b;  // The subsets from the upper half

    size_t m_numSets;     // The number of sets to generate

    size_t m_elementsSum;	// The sum of the elements of m_S
    size_t m_upperBound;  // upper bound for generating sums
    size_t m_lowerBound;  // The current lower bound and


    size_t m_perfectValue; // The value of a perfect partition
    uint64_t m_smallMax;   // The max subset sum of the subsets < perfect
    uint64_t m_largeMin;   // The min subset sum of the subsets > perfect

    bool m_firstCall;     // Double m_numSets after the first call
    const size_t M_MAX_CARDINALITY;	// The max cardinality to consider

    double m_memoryUsage;				// Keep track of memory usage.

  protected :
  	// ---------------------------------------------------------------------------------------
  	// Helper Functions for generating half sets ONLY of cardinality from 0 to MAX_CARDINALITY
  	// TODO: For MAX_CARDINALITY, only include the set if it is in the range of lower and upper
  	// 			 bound since it can only be included in the final HS solution if it is combined
  	//			 with the empty set
  	// ---------------------------------------------------------------------------------------

  	void generateAllSets (const uint64_t S[],
  												vector<SetNodeBitset> &sets,
  												const int first,
  												const int last,
  												const uint64_t cursum,
  												DynamicBitset &curSet,
  												const size_t cardinality);

  	// Helper function, provides missing variables
  	// also sorts sets before returning
  	void generateAllSets(const uint64_t S[], const int N, vector<SetNodeBitset> &sets,
  											 int first, int last, bool sortAscending);


  public :

    HorowitzSahniLowCardinalityMSubsets(const uint64_t S[MAXN], const int N, const int K,
    			const size_t numSets, uint64_t upperBound, const size_t maxCardinality);

    ~HorowitzSahniLowCardinalityMSubsets();
    // smallSets and largeSets are filled in.
    // smallSets contain all of the sets generated that are smaller than perfect
    //    - ordered from largest sum to smallest sum
    // largeSets contain all of the sets generated that are larger than perfect
    //    - ordered from smallest sum to largest sum
    void  generateSetsHS(vector<SetNodeBitset> &smallSets, vector<SetNodeBitset> &largeSets);


    inline double getMemoryUsage() const {
    	return m_memoryUsage;
    }
  };

  typedef HorowitzSahniLowCardinalityMSubsets HSLCMS;
} // end namespace

#endif /* SCHROEPPEL_SHAMIR_CARDINALITY_HPP_ */
