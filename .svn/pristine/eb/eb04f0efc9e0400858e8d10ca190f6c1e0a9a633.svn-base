/*
 * Partition.hpp
 *
 *  Created on: Feb 8, 2014
 *      Author: ethan
 */

#ifndef PARTITION_HPP_
#define PARTITION_HPP_

#include <stdint.h>
#include "../ss/SubsetSum.hpp"	// For DynamicBitset
#include "PartitionUtils.hpp"

// This is an abstract class, uint64_t partition(...) is pure virtual
namespace partition {
class Partition {
protected :
  const vector<uint64_t> m_S;                           // The input elements
  int m_K;                                              // The total number of bins in the problem

  vector<size_t> d_binCounts;                   // For debugging
  size_t m_firstCount;                          // For debugging, how many first sets do we look at?

  Partition(const vector<uint64_t>& S, const int K);
public :
	virtual uint64_t partition(const int K,					// The idx of the current bin we are filling
  								ss::DynamicBitset &elements,		// Bitset with 1 representing elements still left
  								const uint64_t elementsSum,			// The sum of the elements still left
  								const size_t elementsCount,			// The count of the elements left (# of 1's in elements)
  								size_t cardinality,							// The cardinality of sets in the tree we are searching now
  								size_t firstIdx,								// The first idx we can include (in order to stop duplicate permutations)
  								uint64_t partialCost,						// The max sum used so far above this point in the tree
  								uint64_t ub) = 0; 	 						// The maximum sum of a bin determined by the first set sum (k==m_k)


  size_t getFirstCount() const;

  virtual ~Partition();
};

} // end namespace
#endif /* PARTITION_HPP_ */
