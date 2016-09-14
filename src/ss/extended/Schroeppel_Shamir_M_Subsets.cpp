/*
 * Schroeppel_Shamir.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "Schroeppel_Shamir_M_Subsets.hpp"

#include "../../partition/PartitionUtils.hpp"	// for getLowerBound
#include  <iostream>
#include <numeric>  // for accumulate
using std::cout;
using std::endl;
#include "../../utils/MemoryUsage.hpp"

namespace ss {
typedef SchroeppelShamirMSubsets SSMS;

SSMS::SchroeppelShamirMSubsets(const uint64_t S[], const int N, const int K,
                               const size_t numSets, uint64_t upperBound, const uint8_t maxCard) :
  m_K(K),
  m_numSets(numSets),
  m_elementsSum(std::accumulate(S,S+N,(uint64_t) 0)),
  m_lowerBound(0),
  m_upperBound(upperBound),
  m_perfectValue((m_elementsSum + K-1) / K),
  m_smallMax(m_perfectValue-1),
  m_largeMin(m_perfectValue),
  m_firstCall(true),
  m_maxCard(maxCard) {

//  for (int i=0;i<N;i++) {
//    cout << S[i] << endl;
//  }
	m_lowerBound = partition::computeCMin(m_elementsSum, m_K, m_upperBound+1);
  	initializeSets (S, N, m_a, m_b, m_c, m_d, m_maxCard);
//	initializeSets (S, N, m_c, m_d, m_a, m_b, m_maxCard);

}

void  SSMS::generateSetsSS(vector<ss::SetNodeBitset> &smallSets, vector<ss::SetNodeBitset> &largeSets) {


  uint64_t upperBound = m_upperBound;       // Set upper bound to initial upper bound
  uint64_t lowerBound = partition::computeCMin(m_elementsSum, m_K, upperBound);
  if (m_firstCall) {
    m_firstCall = false;      // Do nothing on the first call
  } else {
    m_numSets *= 2;           // For all subsequent, double number of large sets generated
  }

  SSSetMinHeap minHeap;
  SSSetMaxHeap maxHeap;
  deque<SSHeapNode> cdlist;    // list of subset sums removed from CD heap

  SSMinHeap abHeap;    // Min Heap Comparator
  SSMaxHeap cdHeap;    // Max Heap Comparator

  initializeHeaps(m_a, m_b, m_c, m_d, abHeap,cdHeap,upperBound, m_maxCard);
//  initializeHeaps(m_c, m_d, m_a, m_b, abHeap,cdHeap,upperBound, m_maxCard);

  // --------------------------------------------------------------------------
  // Now fill in the allsums vector
  // --------------------------------------------------------------------------

  while (!abHeap.empty()) {                      // until AB heap is empty

    // remove sets with sums above upper bound from cdlist
    while (!cdlist.empty() && abHeap.front().sum() + cdlist.front().sum() > upperBound) {
      cdlist.pop_front();
    }

    // add new sets with sums within bounds to cdlist
    uint64_t topsum = cdHeap.front().sum() + abHeap.front().sum();  //sum of elements on top of heaps

    while (!cdHeap.empty() && topsum >= lowerBound) {       			// new set within bounds

      if (topsum <= upperBound) {                             		// new set within bounds
        cdlist.push_back(cdHeap.front());              							// add heap element to list
      }

      if (cdHeap.front().y() > 0) {                      // there is another element in this column
        SSHeapNode newElement(cdHeap.front().x(),        // keep same row
                              cdHeap.front().y() - 1,    // replace with next combination in column
                              m_c,m_d);                  // sum is c[x] + d[y-1]

        cdHeap.pop();                                    // Remove old head
        cdHeap.push(newElement);                       	 // add new element to heap vector

      } else {                                           // no more elements in this column
        cdHeap.pop();                                  	 // Remove old head
      }

      topsum = cdHeap.front().sum() + abHeap.front().sum();      // sum of two top elements of heaps
    }

    // Combine each element from cdlist with the top element from
    // abHeap, these are all within range
    for (size_t index=0;index<cdlist.size();index++) {

    	uint64_t subsum = cdlist[index].sum() + abHeap.front().sum();  // first subset sum
    	uint8_t card = m_c[cdlist[index] .x()].cardinality() +
    	               m_d[cdlist[index] .y()].cardinality() +
    	               m_a[abHeap.front().x()].cardinality() +
    	               m_b[abHeap.front().y()].cardinality();

    	if (subsum >= lowerBound && subsum <= upperBound &&   // sum is in bounds &&
    	    card <= m_maxCard) {                              // card is not too high
				if (subsum >= m_perfectValue) {

					if ((subsum >= m_largeMin) &&  // Don't generate same subsets generated on previous runs
					    (maxHeap.size() < m_numSets || subsum < maxHeap.front().sum())) {  // Don't have enough yet or smaller than largest so far

						maxHeap.push(SetNodeBitset(subsum, m_c[cdlist[index].x()].set() |
																					     m_d[cdlist[index].y()].set() |
																					     m_a[abHeap.front().x()].set() |
																					     m_b[abHeap.front().y()].set()));
					}

					// What about duplicates?
					while (maxHeap.size() > m_numSets) {
						maxHeap.pop();

						upperBound = maxHeap.front().sum();
						lowerBound = partition::computeCMin(m_elementsSum, m_K, upperBound);

						while (!minHeap.empty() && minHeap.front().sum() < lowerBound) {
							minHeap.pop();
						}
					}
				} else if (subsum <= m_smallMax) { // Don't generate same subsets generated on previous runs

				  minHeap.push(SetNodeBitset(subsum, m_c[cdlist[index].x()].set() |
																				     m_d[cdlist[index].y()].set() |
																				     m_a[abHeap.front().x()].set() |
																				     m_b[abHeap.front().y()].set()));
				}
    	}
    }

    // If there is another element from the b set, add it.
    // Then, remove top element from abHeap.
    if (abHeap.front().y() < m_b.size() - 1) {          // there is another element in this column
      SSHeapNode newNode(abHeap.front().x(),            // new element being added to heap
                       abHeap.front().y() + 1,        // replace with next combination in this column
                       m_a,m_b);

      abHeap.pop();                            // Pop old value
      abHeap.push(newNode);                    // add new element to heap vector
    } else {
      abHeap.pop();
    }
  }

 // process_mem_usage(m_memoryUsage);
  //cout << "SS Memory: " << m_memoryUsage << endl;

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
//  m_upperBound = (2 * (m_upperBound - m_perfectValue)) + m_perfectValue;    // This guesses the next upper/lower bound
//  m_lowerBound = partition::computeCMin(m_elementsSum, m_K, m_upperBound+1);  // Need to get rid of local lowerBound and upperBound variables to use

  size_t numSmallSets = minHeap.size();
  size_t numLargeSets = maxHeap.size();
	
	if (numLargeSets == 0) {
		exit(0);
	}
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

//  { // TMP DEBUG for Cardinality HS
//		cout << "Lower Bound: " << lowerBound << endl
//				 << "Upper Bound: " << upperBound << endl << endl;
//
//		SimpleTimer timer;
//		cout << "Running HS...";
//		cout.flush();
//		vector<ss::SSSetNode> allSets;
//		ss::generateSetsHSCardinality(&m_S[0],m_S.size(),lowerBound,upperBound,6,allSets);
//		cout << "Done" << endl
//				 << "Time : " << timer.timeElapsed() << " s" << endl
//				 << "Count: " << allSets.size() << endl << endl;
//  }

//  // DEBUG START
//  cout << "Num Large  : " << largeSets.size() << endl
//  		 << "Num Small  : " << smallSets.size() << endl
//  		 << "Lower Bound: " << lowerBound << endl
//  		 << "Upper Bound: " << upperBound << endl
//  		 << "Num Sets   : " << m_numSets << endl;
//  uint64_t *tmp = new uint64_t[m_S.size()];
//  for (size_t i=0;i<m_S.size();i++) {
//  	tmp[i] = m_S[i];
//  }
//
//  for (size_t i=0;i<largeSets.size();i++) {
////  	if (largeSets[i].set().test(0)) {
//    if (i % 100 == 0) {
//  		cout << largeSets[i].toString(tmp) << endl;
//  	}
//  }

//  cout << "Count: " << smallSets.size() + largeSets.size() << endl;
//  for (size_t i=0;i<largeSets.size();i+=1000) {
//  	cout << largeSets[i].toString(tmp) << endl;
//  }
//  cout << largeSets.back().toString(tmp) << endl;
//  delete []tmp;
//exit(0);
// DEBUG END
}


} // end namespace
