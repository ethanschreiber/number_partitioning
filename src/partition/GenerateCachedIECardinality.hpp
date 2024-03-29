/*
 * CachedIECard.hpp
 *
 *  Created on: Feb 7, 2014
 *      Author: ethan
 */

#ifndef CACHEDIECARDINALITY_HPP_
#define CACHEDIECARDINALITY_HPP_
#include "Partition.hpp"
#include "CachedIETrees.hpp"
#include "../Utils.hpp"
// ============================================================================
//
// C a c h e d I E C a r d i n a l i t y
//
// ============================================================================
namespace partition {

class GenerateCachedIECardinality  {
private :
	bool m_isTmpVerbose;
  CachedIETrees *m_cachedIETrees;   // The cached IE data structure;
  Partition *m_partition; 					// Pointer to Partition instance
  const vector<uint64_t> m_S; 			// The input elements
  int m_K;													// The number of subsets in a partition

  // Debugging variables

  vector<size_t> d_timesCalled;			// The # of times generation is called for a particular cardinality
  vector<size_t> d_nodesExplored;		// The # of nodes generated over all generation calls
protected :

  uint64_t generateRestCIE(const CachedIENode *node,	  // The current node in the IE tree
  								         const int K,								  // The idx of the current bin we are filling
  								         ss::DynamicBitset &elements, // Bitset with 1 representing elements still left
  								         const uint64_t elementsSum,	// The sum of the elements still left
  								         const size_t elementsCount,	// The count of the elements left (# of 1's in elements)
  								         const uint64_t pathSum,      // The sum of the elements inclued so far in this tree recursion
  								         size_t cardinality,					// The cardinality of sets in the tree we are searching now
  								         size_t firstIdx,							// The first idx we can include (in order to stop duplicate permutations)
  								         size_t firstIdxMax,					// The last idx we can include so we can still make enough completions with this cardinality
  								         uint64_t lb,								  // The minimum sum of a set that could lead to a solution
  								         uint64_t partialCost,				// The max sum used so far above this point in the tree
  								         uint64_t ub); 					      // The maximum sum of a bin determined by the first set sum (k==m_k)


  uint64_t generateRestIE(const size_t elementIdx,      // The idx of the next element to include/exclude
                          const int K,                  // The idx of the current bin we are filling
                          ss::DynamicBitset &elements,  // Bitset with 1 representing elements still left
                          const uint64_t elementsSum,   // The sum of the elements still left
                          const size_t elementsCount,   // The count of the elements left (# of 1's in elements)
                          const uint64_t pathSum,       // The sum of the included elements on the current path
                          const size_t pathCardinality, // The number of included elements on the current path
                          const size_t cardinality,     // The cardinality of sets in the tree we are searching now
                          size_t firstIdx,              // The first idx we can include (in order to stop duplicate permutations)
                          const size_t firstIdxMax,     // The last idx we can include so we can still make enough completions with this cardinality
                          const uint64_t lb,            // The minimum sum of a set that could lead to a solution
                          const uint64_t partialCost,   // The max sum used so far above this point in the tree
                          uint64_t ub);                 // The maximum sum of a bin determined by the first set sum (k==m_k)


  bool canPrune(size_t elementIdx,      // The idx of the next element to include/exclude
  							const ss::DynamicBitset &elements, // Bitset with 1 representing elements still left
  							const uint64_t pathSum,      // The sum of the elements inclued so far in this tree recursion
      					const size_t pathCardinality, // The number of included elements on the current path
      					const size_t cardinality,     // The cardinality of sets in the tree we are searching now
      					const uint64_t lb,            // The minimum sum of a set that could lead to a solution
      					const uint64_t ub);                 // The maximum sum of a bin determined by the first set sum (k==m_k)


  void push_back(const ss::SetNodeBitset &node);

public :
  GenerateCachedIECardinality(Partition *partition, const vector<uint64_t> &S, CachedIETrees *cachedIETrees,
  													uint64_t elementsSum, int K, uint64_t upperBound, size_t numSet);
  ~GenerateCachedIECardinality();

  void setTmpVerbose() {
	  m_isTmpVerbose = true;
  }

  uint64_t generate(const int K,									// The idx of the current bin we are filling
  								 ss::DynamicBitset &elements,	// Bitset with 1 representing elements still left
  								 const uint64_t elementsSum,	// The sum of the elements still left
  								 const size_t elementsCount,	// The count of the elements left (# of 1's in elements)
  								 size_t cardinality,					// The cardinality of sets in the tree we are searching now
  								 size_t firstIdx,							// The first idx we can include (in order to stop duplicate permutations)
  								 uint64_t lb,									// The minimum sum of a set that could lead to a solution
  								 uint64_t partialCost,				// The max sum used so far above this point in the tree
  								 uint64_t ub); 								// The maximum sum of a bin determined by the first set sum (k==m_k)

	string cardinalityCountsToString() const;
	void resetCardinalityCounts() {
    for (size_t i=0;i<m_S.size();i++) {
      d_timesCalled[i] = 0;
      d_nodesExplored[i]=0;
    }
	}
};

typedef GenerateCachedIECardinality GCIEC;
}

#endif /* CACHEDIECARD_HPP_ */
