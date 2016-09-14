/*
 * Schroeppel_Shamir.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "HorowitzSahniLowCardinalityMSubsets.hpp"

#include "../../partition/PartitionUtils.hpp"	// for getLowerBound
#include "../SchroeppelShamir.hpp"
#include "../../utils/MemoryUsage.hpp"

#include  <iostream>
#include <numeric>  // for accumulate
using std::cout;
using std::endl;

namespace ss {

HSLCMS::HorowitzSahniLowCardinalityMSubsets(const uint64_t S[], const int N, const int K,
																						const size_t numSets, uint64_t upperBound, const size_t maxCardinality) :
  m_S(S,S+N), m_K(K), m_numSets(numSets), m_elementsSum(std::accumulate(S,S+N,(uint64_t) 0)),
  m_upperBound(upperBound), m_lowerBound(partition::computeCMin(m_elementsSum, m_K, m_upperBound+1)),
  m_perfectValue((m_elementsSum + K-1) / K), m_smallMax(m_perfectValue-1),m_largeMin(m_perfectValue),
  m_firstCall(true), M_MAX_CARDINALITY(maxCardinality){

//	cout << "Generate A:" << endl;
	generateAllSets(S, N, m_a, 0    , N/2-1  , true);   // combinations of 1st half of numbers (sorted ascending)

//	cout << "\nGenerate B:" << endl;
	generateAllSets(S, N, m_b, N/2  , N-1    , false);  // combinations of 2nd quarter of numbers (sorted descending)



//	cout << "Upper Bound: " << upperBound << endl;
//	cout << "A: ";
//	for (size_t i=0;i<m_a.size();i++) {
//		cout << m_a[i].toString(&m_S[0]) << endl;
//	}
//
//	cout << endl;
//
//	cout << "B: ";
//	for (size_t i=0;i<m_b.size();i++) {
//		cout << m_b[i].toString(&m_S[0]) << endl;
//	}


//	cout << "A size: " << m_a.size() << endl
//			 << "B size: " << m_b.size() << endl;
}

ss::HorowitzSahniLowCardinalityMSubsets::~HorowitzSahniLowCardinalityMSubsets() {

}

void  HSLCMS::generateSetsHS(vector<ss::SetNodeBitset> &smallSets, vector<ss::SetNodeBitset> &largeSets) {

  SSSetMinHeap minHeap;
  SSSetMaxHeap maxHeap;

  uint64_t upperBound = m_upperBound;       // Set upper bound to initial upper bound
  uint64_t lowerBound = partition::computeCMin(m_elementsSum, m_K, upperBound);

	size_t bAnchor = 0;                              			// bPtr starts searching from here

	for (size_t aPtr=0; aPtr < m_a.size(); aPtr++) {
		const uint64_t &aValue = m_a[aPtr].sum();        			// For readability, the sum of the a subset.
		const size_t aCardinality = m_a[aPtr].set().count();	// The cardinality of the a subset.

		while (bAnchor < m_b.size() &&                   			// Move bAnchor to highest value that might fit
					 aValue + m_b[bAnchor].sum() > upperBound) {    // Keep looking at smaller number in b until we
			bAnchor++;                                   				// are under upper bound. (remember b is descending)
		}

		// Now add all aValue + b[bPtr] that are within range
		for (size_t bPtr=bAnchor;(bPtr < m_b.size()) && (aValue + m_b[bPtr].sum() >= lowerBound); bPtr++) {

			// Only add if subset cardinality is no greater than MAX_CARDINALITY
			if (aCardinality + m_b[bPtr].set().count() <= M_MAX_CARDINALITY) {
				uint64_t subsum = aValue + m_b[bPtr].sum();

				if (subsum >= m_perfectValue) {

					if ((subsum >= m_largeMin) &&  // Don't generate same subsets generated on previous runs
					    (maxHeap.size() < m_numSets || subsum < maxHeap.front().sum())) {  // Don't have enough yet or smaller than largest so far
						maxHeap.push(SetNodeBitset(subsum,
																				m_a[aPtr].set() | m_b[bPtr].set()));
					}

					// TODO: What about duplicates?
					while (maxHeap.size() > m_numSets) {
						maxHeap.pop();

						upperBound = maxHeap.front().sum();
						lowerBound = partition::computeCMin(m_elementsSum, m_K, upperBound);

						while (!minHeap.empty() && minHeap.front().sum() < lowerBound) {
							minHeap.pop();
						}

					}
				} else if (subsum <= m_smallMax){

					minHeap.push(SetNodeBitset(subsum,
																			m_a[aPtr].set() | m_b[bPtr].set()));
				}
			}
		}
	}

	process_mem_usage(m_memoryUsage);
//	cout << "HS Low Card Memory: " << m_memoryUsage << endl;
  // If we run this function again, we want to generate sets with sum
  // < m_smallMax and > m_largeMin; On a number line:
  //
  // So for each call, we generate values between lowerBound and m_smallMax
  // as well as between largeMin and upperBound
  // lowerBound      m_smallMax     m_perfectValue         m_largeMin        upperBound
  //     |------------------|------------------|------------------|------------------|

  if (!minHeap.empty()) {
    m_smallMax = minHeap.front().sum() - 1;
  }
  if (!maxHeap.empty()) {
    m_largeMin = maxHeap.front().sum() + 1;
  }

  size_t numSmallSets = minHeap.size();
  size_t numLargeSets = maxHeap.size();

  smallSets.resize(numSmallSets);
  largeSets.resize(numLargeSets);

  for (size_t i=0;i<numSmallSets;i++) {
    smallSets[numSmallSets-i-1] = minHeap.front();
    minHeap.pop();
  }

  for (size_t i=0;i<numLargeSets;i++) {
    largeSets[numLargeSets-i-1] = maxHeap.front();
    maxHeap.pop();
  }
}

// ---------------------------------------------------------------------------------------
// Helper Functions for generating half sets ONLY of cardinality from 0 to MAX_CARDINALITY
// TODO: For MAX_CARDINALITY, only include the set if it is in the range of lower and upper
// 			 bound since it can only be included in the final HS solution if it is combined
//			 with the empty set
// ---------------------------------------------------------------------------------------

void HSLCMS::generateAllSets (const uint64_t S[],
															vector<SetNodeBitset> &sets,
															const int first,
															const int last,
															const uint64_t cursum,
															DynamicBitset &curSet,
															const size_t cardinality)

{
	if (first > last) {                                    // set is completed
		if (cardinality == M_MAX_CARDINALITY) {
			if (cursum >= m_lowerBound && cursum <= m_upperBound) {		// If at max, must be within bounds
				sets.push_back(SetNodeBitset(cursum,curSet));
			}
		} else {

			if (cursum <= m_upperBound) {														// If less than max, has to be less than or equal to upper
				sets.push_back(SetNodeBitset(cursum,curSet));
			}
		}
	} else {                                               // set not yet completed

		// exclude next element
		generateAllSets (S, sets, first+1, last, cursum, curSet, cardinality);


		if (cardinality < M_MAX_CARDINALITY) {	// If less than MAX_CARDINALITY
			// Include next element
			curSet[first] = true;
			generateAllSets (S, sets, first+1, last, cursum + S[first], curSet, cardinality+1);
			curSet[first] = false;
		}
	}
}

// Helper function, provides missing variables
// also sorts sets before returning
void HSLCMS::generateAllSets(const uint64_t S[], const int N, vector<SetNodeBitset> &sets,
										 	 	 	 	 int first, int last, bool sortAscending) {
	DynamicBitset curSet(N);                                   // The number of total input elements
	generateAllSets(S,sets,first,last,0ll,curSet,0);

	if (sortAscending) {
		std::sort( sets.begin(), sets.end() ,std::less<SetNodeBitset>());  // Sort in ascending order
	} else {
		std::sort( sets.begin(), sets.end() ,std::greater<SetNodeBitset>());  // Sort in descending order
	}

}

} // end namespace

