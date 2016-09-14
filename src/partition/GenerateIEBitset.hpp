/*
 * MoffittBitSet.hpp
 *
 *  Created on: Aug 20, 2013
 *      Author: ethan
 */

#ifndef MOFFITT_BIT_SET_HPP_
#define MOFFITT_BIT_SET_HPP_

#include "Partition.hpp"

namespace partition {

// ============================================================================
//
// G e n e r a t e I E B i t s e t - This is simple exclusion/exclusion which
// uses a bitset to keep track of the remaining elements in S.
//
// ============================================================================
class GenerateIEBitset {
private :
  Partition *m_partition; 		// Pointer to Partition instance
  const vector<uint64_t> m_S; // The input elements
  int m_K;                    // The total number of bins in the problem

protected :
  uint64_t generateRest(const int K,									// The idx of the current bin we are filling
    								 	 	ss::DynamicBitset &elements,	// Bitset with 1 representing elements still left
    								 	 	const uint64_t elementsSum,		// The sum of the elements still left
												const size_t elementsCount,		// The count of the elements left (# of 1's in elements)
                        const uint64_t includedSum, 	// The sum of the included elements on the current path
                        const size_t elementIdx,    	// The idx of the next element to include/exclude
                        size_t cardinality,						// Just for passing back to Partition
    								 	 	uint64_t lb,									// The minimum sum of a set that could lead to a solution
    								 	 	uint64_t partialCost,					// The max sum used so far above this point in the tree
    								 	 	uint64_t ub); 								// The maximum sum of a bin determined by the first set sum (k==m_k)

public :

  GenerateIEBitset(Partition *partition,
									const vector<uint64_t> &S,
									int K);

  uint64_t generate(const int K,									// The idx of the current bin we are filling
    								 ss::DynamicBitset &elements,	// Bitset with 1 representing elements still left
    								 const uint64_t elementsSum,	// The sum of the elements still left
    								 const size_t elementsCount,	// The count of the elements left (# of 1's in elements)
    								 size_t cardinality,					// Just for passing back to Partition
    								 uint64_t lb,									// The minimum sum of a set that could lead to a solution
    								 uint64_t partialCost,				// The max sum used so far above this point in the tree
    								 uint64_t ub); 								// The maximum sum of a bin determined by the first set sum (k==m_k)
};

} // End Namespace

#endif /* MOFFITT_HPP_ */
