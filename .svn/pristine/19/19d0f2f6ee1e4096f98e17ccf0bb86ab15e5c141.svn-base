/*
 * PartitionIterativeWeakeningLowCardinality.cpp
 *
 *  Created on: Jul 31, 2013
 *      Author: ethan
 */

#include "PartitionIterativeWeakeningLowCardinality.hpp"
#include "../ss/extended/Extended_Schroeppel_Shamir.hpp"
#include "../utils/MemoryUsage.hpp"

namespace partition {

// ===================================================================================
//
// C a c h e d I t e r a t i v e W e a k e n i n g L o w C a r d i n a l i t y (PIWLC)
//
// ====================================================================================

PIWLC::PartitionIterativeWeakeningLowCardinality(const vector<uint64_t>& S,
    uint64_t elementsSum, int K, uint64_t ub, size_t numSets) :
    Partition(S, K), m_maxCardinality(((S.size() + K - 1) / K)+1), m_cachedIETrees(
        new CachedIETreesLowCardinalityHS(S, K, ub, numSets, m_maxCardinality)), m_cachedIE(
        new GenerateCachedIECardinality(this, S, m_cachedIETrees, elementsSum,
            K, ub, numSets)) {

}

PIWLC::~PartitionIterativeWeakeningLowCardinality() {
  delete m_cachedIE;
  delete m_cachedIETrees;
}

// ----------------
// Search Functions
// ----------------

// Function used to complete all subsequent bins for CIE trees
uint64_t PIWLC::partition(const int K, // The idx of the current bin we are filling
    ss::DynamicBitset &elements, // Bitset with 1 representing elements still left
    const uint64_t elementsSum,	// The sum of the elements still left
    const size_t elementsCount,	// The count of the elements left (# of 1's in elements)
    size_t cardinality,	// The cardinality of sets in the tree we are searching now
    size_t firstIdx,// The first idx we can include (in order to stop duplicate permutations)
    uint64_t partialCost,	// The max sum used so far above this point in the tree
    uint64_t ub) { 				// The best solution found so far

  d_binCounts[K]++;
  uint64_t currentValue;
  if (K == 1) {                          			// If we have reached the last bin
    currentValue = elementsSum;							// All elements must go in it
  } else if (elementsCount == 0) { // If there are no elements left (this should never happen)
    currentValue = 0;	// Then all remaining bins are empty and they have to fit
  } else {
    uint64_t lb = computeCMin(elementsSum, K, ub);// Minimum value for remaining elements

    if (lb < ub) {

//    	if (K >= 9) {
//    		cout << "CV Before: " << currentValue  << " K: " << K << endl;
//    	}
      currentValue = m_cachedIE->generate(K, elements, elementsSum,
          elementsCount, cardinality, firstIdx, lb, partialCost, ub);
//  	if (K >= 9) {
//  		cout << "CV After : " << currentValue  << " K: " << K << endl;
//  		cout << "CV Partia: " << partialCost  << endl << endl;
//  	}
    } else {
      currentValue = ub;
    }
  }

  return std::max(partialCost, currentValue);
}

uint64_t PIWLC::partitionFirst(const int N, const uint64_t elementsSum, uint64_t upperBound) {
  uint64_t currentValue = 1; 							// Keep going until we find a solution

  uint64_t firstSum = 0;									// With firstSum being the largest sum

  SimpleTimer timer;
  while (currentValue > firstSum) {	// Keep going until we have found an optimal solution

    m_cachedIETrees->addSetsToCache();    // Add necessary values to prefix tree

    firstSum = (m_cachedIETrees->getMax())->sum();	// The sum of the first set


    const ss::SetNodeBitset &setNode =	// The next set to consider for the first bin
        *(m_cachedIETrees->getMax());
    cout << "First Sum          : " << firstSum << "   Card: " << (int) setNode.cardinality() << endl;
    if (setNode.set().count() <= m_maxCardinality) {// This is arbitrary, should be <= but probably shouldn't exist

//      cout << "First Sum          : " << std::setw(3) << m_firstCount << ": " << firstSum << endl;
      m_firstCount++;
      ss::DynamicBitset elements = ~setNode.set();// All that aren't used are left


      m_cachedIE->resetCardinalityCounts();

      currentValue = partition(m_K - 1, elements, elementsSum - firstSum,
          elements.count(), m_cachedIETrees->getMinCardinality(), 0, firstSum,
          firstSum+1);
//      if (currentValue < upperBound) {
//        upperBound = currentValue;
//        cout << "New Upper Bound    : " << upperBound << endl;
//      }
//      cout << m_cachedIE->cardinalityCountsToString();
//      cout << "----------------------------------------------------------" << endl;
    } else {
      currentValue = firstSum + 1;
    }

  }

// FINAL CHECK, NEED THIS FOR OPTIMALITY
  timer.reset();
  ss::DynamicBitset elements(N);
  for (int i = 0; i < N; i++) {
    elements[i] = true;
  }
  m_cachedIE->resetCardinalityCounts();

  currentValue = partition(m_K,			        // The idx of the current bin we are filling
      elements,								              // Bitset with 1 representing elements still left
      elementsSum,								          // The sum of the elements still left
      elements.count(),	                    // The count of the elements left (# of 1's in elements)
      m_cachedIETrees->getMinCardinality(), // The cardinality of sets in the tree we are searching now
      0,                                    // The first idx we can include (in order to stop duplicate permutations)
      0, 								                    // The max sum used so far above this point in the tree
      currentValue);							          // The best solution found so far

  return currentValue;
}

size_t PIWLC::getSsCalls() const {
  return m_cachedIETrees->getSsCalls();
}

double PIWLC::getSsTime() const {
  return m_cachedIETrees->getSsTime();
}
// ===============================
// E x e c u t e   F u n c t i o n
// ===============================

uint64_t executeCIWLowCardinality(const PartitionProblem &problem,
                                  ProblemStats &stats, size_t numSets) {


  uint64_t upperBound = getUpperBound(problem);

  // SimpleTimer timer;
  PartitionIterativeWeakeningLowCardinality mc(
      vector<uint64_t>(problem.S, problem.S + problem.N), problem.sum,
      problem.K, upperBound, numSets);

//  cout << "Construction Time: " << timer.timeElapsed() << endl;
  uint64_t returnValue = mc.partitionFirst(problem.N, problem.sum, upperBound);

  stats.firstCount = mc.getFirstCount();
  stats.ssTime = mc.getSsTime();
  stats.ssCalls = mc.getSsCalls();
  process_mem_usage(stats.residentMemory);
  return returnValue;
}

} // End Namespace
