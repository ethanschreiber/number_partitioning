/*
 * MoffittCachedCardinality.hpp
 *
 *  Created on: Aug 20, 2013
 *      Author: ethan
 */

#ifndef CACHED_ITERATIVE_WEAKENING_HPP_
#define CACHED_ITERATIVE_WEAKENING_HPP_


#include "GenerateCachedIECardinality.hpp"

namespace partition {

// ====================================
// The main function to call externally
// ====================================

uint64_t executeCIW(const PartitionProblem &problem, ProblemStats &stats,size_t numSets);

// ============================================================================
//
// P a r t i t i o n I t e r a t i v e W e a k e n i n g
//
// ============================================================================

class PartitionIterativeWeakening : public Partition {
protected :
	CachedIETreesAllCardinalitySS *m_cachedIETrees;	// The Cached IE Trees class
	GenerateCachedIECardinality *m_cachedIE;    		// The Cached IE class (For searching m_cachedIETrees)
	size_t m_numSetsToAdd;													// Add 1 set then 2 sets then 4 sets etc.
	uint64_t m_upperBound; // Keep track of upper bound in case optimal == approximate solution	
public :

  PartitionIterativeWeakening(const vector<uint64_t> &S,	// The input elements
  												 uint64_t elementsSum,					// The sum of the input elements
  												 int K,													// The number of bins to partition S into
  												 uint64_t upperBound, 					// A known upper bound for the optimal bin capacity
  												 size_t numSets);								// The # of sets to generate for each SS call

  virtual ~PartitionIterativeWeakening();

  // Function used to complete the first bin
  uint64_t partitionFirst(const int N,										// The number of input elements
  												const uint64_t elementsSum);		// The sum of the input elements

    // Function used to complete all subsequent bins
  uint64_t partition(const int K,													// The idx of the current bin we are filling
  									 ss::DynamicBitset &elements,					// Bitset with 1 representing elements still left
  									 const uint64_t elementsSum,					// The sum of the elements still left
  									 const size_t elementsCount,					// The count of the elements left (# of 1's in elements)
  									 size_t cardinality,									// The cardinality of sets in the tree we are searching now
  									 size_t firstIdx,											// The first idx we can include (in order to stop duplicate permutations)
  									 uint64_t maxSoFar,										// The max sum used so far above this point in the tree
  									 uint64_t bestSoFar); 								// The maximum sum of a bin determined by the first set sum (k==m_k)

  // Debugging functions
  size_t getSsCalls() const;
	double getSsTime() const;

};

typedef PartitionIterativeWeakening PIW; // Full name too long
} // End Namespace

#endif /* MOFFITT_HPP_ */
