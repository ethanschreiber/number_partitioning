/*
 * PartitionIterativeWeakening.cpp
 *
 *  Created on: Jul 31, 2013
 *      Author: ethan
 */

#include "PartitionIterativeWeakening.hpp"
#include "../utils/MemoryUsage.hpp"
#include <cstring>
#include <algorithm>
namespace partition {



// ============================================================================
//
// P a r t i t i o n I t e r a t i v e W e a k e n i n g (PIW)
//
// ============================================================================
PIW::PartitionIterativeWeakening(const vector<uint64_t>& S, uint64_t elementsSum, int K, uint64_t upperBound, size_t numSets)
: Partition(S,K),
  m_cachedIETrees(new CachedIETreesAllCardinalitySS(S,K,upperBound,numSets)),
  m_cachedIE(new GenerateCachedIECardinality(this,S,m_cachedIETrees,elementsSum,K,upperBound, numSets))
, m_numSetsToAdd(1), m_upperBound(upperBound)
{
}

PIW::~PartitionIterativeWeakening() {
	delete m_cachedIETrees;
	delete m_cachedIE;
}
// ----------------
// Search Functions
// ----------------
uint64_t PIW::partition(const int K,								// The idx of the current bin we are filling
										 ss::DynamicBitset &elements,		// Bitset with 1 representing elements still left
										 const uint64_t elementsSum,		// The sum of the elements still left
										 const size_t elementsCount,		// The count of the elements left (# of 1's in elements)
										 size_t cardinality,						// The cardinality of sets in the tree we are searching now
										 size_t firstIdx,								// The first idx we can include (in order to stop duplicate permutations)
										 uint64_t partialCost,					// The max sum used so far above this point in the tree
										 uint64_t ub) { 								// The best solution found so far

	d_binCounts[K]++;
	uint64_t currentValue;

//		cout << "------------------------------------------------------------------" << endl;
//		cout << "Partition K: " << K << " MAX K: " << m_K << "  EL: " << elementsCount << "  CARD: " << cardinality << endl;

  if (K==1) {                          							// If we have reached the last bin
  	currentValue = elementsSum;							  			// All elements must go in it
  } else if (!elements.any()) {        							// If there are no elements left (this should never happen)
    currentValue = 0;			                					// Then all remaining bins are empty and they have to fit
  } else {
    uint64_t lb = computeCMin(elementsSum,K,ub);		// Minimum value for remaining elements
    if (lb < ub) {
      currentValue = m_cachedIE->generate(K,elements,elementsSum, elementsCount,
      																	cardinality, firstIdx, lb, partialCost, ub);
    } else {
    	currentValue = ub;
    }
  }


  return std::max(partialCost,currentValue);
}


uint64_t PIW::partitionFirst(const int N,
                        const uint64_t elementsSum) {
  uint64_t currentValue=1;												// Keep going until we find a solution
  uint64_t firstSum = 0;													// With firstSum being the largest sum

//	int tmpCount = 1;
  while (currentValue > firstSum) {							  // Keep going until we have found an optimal solution

  	m_cachedIE->resetCardinalityCounts();
  	SimpleTimer tmpTimer;


//  	// First add one set above perfect, then 2 sets, then 4 sets, then 8 sets...
//  	// This avoids a huge number of failures in a row but it means that
//  	// when we do find a solution with the largest included, it is not necessarily
//  	// optimal, so one last run must be done at the end
//		bool hasMore=true;
//  	for (int i=0;i<m_numSetsToAdd && hasMore;i++) {
//  		hasMore = m_cachedIETrees->addSetsToCache();          	// Add necessary values to prefix tree
//  	}m_numSetsToAdd *= 2;

  	bool hasMore = m_cachedIETrees->addSetsToCache();          	// Add necessary values to prefix tree
  	if (!hasMore) {
  		cout << "HAVE NO MORE, EXITING" << endl;
  		exit(0);
  	}

  	const ss::SetNodeBitset &setNode =								// The next set to consider for the first bin
    		*(m_cachedIETrees	->getMax());
    firstSum = setNode.sum();  										// The sum of this set
		
		if (firstSum == m_upperBound) {
			return firstSum;
		}
//		cout << "First Sum          : " << firstSum << "   Card: " << (int) setNode.cardinality() << endl;
					
    m_firstCount++;
    ss::DynamicBitset elements = ~setNode.set();	// All that aren't used are left
		
    currentValue = partition(m_K-1,elements,elementsSum-firstSum, elements.count(),
														 m_cachedIETrees->getMinCardinality(), 0, firstSum, firstSum+1);  // Perform search

//    cout << "Current Value: " << currentValue << endl;
//    if (!hasMore) {
//      break;
//    }
//    if (setNode.set().count() <= 5) {
//		 cout << std::setw(2) << tmpCount << ": " << tmpTimer.timeElapsed();
//		 tmpCount++;
//		 cout << m_cachedIE->cardinalityCountsToString();
//		 cout << "------------------------------------------------" << endl << endl;
//    }

#if PAPER_DEBUG == 1
		const string CARD_TREES_FILENAME("/home/ethan/workspace/aaai14/card_trees.dot");
		cout << "---------------------" << endl
		     << "FC: " << getFirstCount() << endl
		     << "Writing " << CARD_TREES_FILENAME << endl
		     << "---------------------" << endl;
		std::ofstream cardTreesFile(CARD_TREES_FILENAME);
		cardTreesFile << m_cachedIE->toDot() << endl;
		cardTreesFile.close();
#endif

  }

#if PAPER_DEBUG == 1
  uint64_t CMin = computeCMin(elementsSum,m_K,currentValue+1); // Minimum value for remaining elements
  cout << "Sets: " << endl;

  const string SETS_FILENAME("/home/ethan/workspace/aaai14/sets.tex");
  cout << "Writing Sets to " << SETS_FILENAME << endl;
  std::ofstream setsFile(SETS_FILENAME);
  setsFile << m_cachedIE->printSets(CMin,currentValue) << endl;
  setsFile.close();
#endif


  // FINAL CHECK, NEED THIS FOR OPTIMALITY
//  if (m_numSetsToAdd > 1) {														// If we added more than one set at a time, we need to double check for optimality
////  	cout << "Num Sets To Add: " << m_numSetsToAdd << endl;
//			SimpleTimer timer;
//			ss::DynamicBitset elements(N);
//			for (size_t i=0;i<N;i++) {
//				elements[i] = true;
//			}
//			currentValue = partition(m_K,elements,elementsSum, elements.count(),
//															 m_cachedIETrees->getMinCardinality(), 0, 0, currentValue);
////			cout << "Final: " << timer.timeElapsed() << endl;
//  }


    return currentValue;

}


size_t PIW::getSsCalls() const {
	return m_cachedIETrees->getSsCalls();
}

double PIW::getSsTime() const {
	return m_cachedIETrees->getSsTime();
}

// ===============================
// E x e c u t e   F u n c t i o n
// ===============================


uint64_t executeCIW(const PartitionProblem &problem, ProblemStats &stats, size_t numSets) {

  uint64_t upperBound = getUpperBound(problem);
//    SimpleTimer timer;
	uint64_t *SCopy = new uint64_t[problem.N];
	memcpy(SCopy,problem.S,sizeof(uint64_t) * problem.N);

//	cout << "# Sets: " << numSets << endl;
	std::sort( SCopy, SCopy+problem.N ,std::greater<uint64_t>());

  PartitionIterativeWeakening piw(vector<uint64_t>(SCopy,SCopy+problem.N),
                          problem.sum, problem.K,upperBound, numSets);

//  cout << "Construction Time: " << timer.timeElapsed() << endl;
  uint64_t returnValue = piw.partitionFirst(problem.N, problem.sum);

  stats.firstCount = piw.getFirstCount();

  stats.ssTime = piw.getSsTime();
  stats.ssCalls = piw.getSsCalls();

  process_mem_usage(stats.residentMemory);
  delete [] SCopy;
  return returnValue;
}

} // End Namespace
