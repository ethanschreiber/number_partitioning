/*
 * Schroeppel_Shamir.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "SchroeppelShamir.hpp"
#include "HorowitzSahni.hpp"
#include "KK.hpp"
#include "../utils/MemoryUsage.hpp"
#include "../Utils.hpp"
#include  <iostream>
#include <iomanip>
using std::cout;
using std::endl;

namespace ss {


// ============================================================================
//
// Schroeppel and Shamir Helper Functions *** SETS ***
//
// ============================================================================

void initializeSets(const uint64_t S[MAXN], const int n,
                    vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                    const uint8_t maxCardinality) {

  generateAllSets(S, n, a, 0    , n/4-1  ,true,maxCardinality); // combinations of 1st quarter of numbers
  generateAllSets(S, n, b, n/4  , n/2-1  ,true,maxCardinality); // combinations of 2nd quarter of numbers
  generateAllSets(S, n, c, n/2  , n*3/4-1,true,maxCardinality); // combinations of 3rd quarter of numbers
  generateAllSets(S, n, d, n*3/4, n-1    ,true,maxCardinality); // combinations of 4th quarter of numbers
}

void initializeHeaps(vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                     SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper, const uint8_t maxCardinality) {

  // --------------------------------------------------------------------------
  //construct min heap of initial combinations of sets from array A and array B
  // --------------------------------------------------------------------------

  // The crdinalities always work here since the cardinalities from a
  // are all <= maxCardinality and we take the empty set from b
  for (size_t aIndex = 0; aIndex < a.size() && a[aIndex].sum() <= upper; aIndex++) {  // for each set in A array
    abheap.vector_push_back(SSHeapNode(aIndex,      // aptr pointer is same as index
                                0,           // smallest value of bsums is first element
                                a[aIndex])); // b.front().sum == 0 --> [a[index].sum + b.front().sum] == a[index].sum,
  }

  abheap.make_heap();


  // --------------------------------------------------------------------------
  //construct max heap of initial combinations of sets from array C and array D
  // --------------------------------------------------------------------------

  int dIndex = d.size() - 1;                  // index to largest value in dsums

  for (size_t cIndex = 0; cIndex < c.size(); cIndex++) {  // process csums in increasing order
    if (c[cIndex].sum() > upper) {
      break;                                    // no more entries to add to heap
    }

    // Search through d until we find the first set whose sum is not too large
    while (c[cIndex].sum() + d[dIndex].sum() > upper) {
      dIndex--;
    }

    // Now make sure the cardinalities work too, but don't change dIndex,
    // change dIndexCard since we don't want this to affect the next c[cIndex]
    int dIndexCard = dIndex;
    while (c[cIndex].cardinality() + d[dIndexCard].cardinality() > maxCardinality) {
      dIndexCard--;
    }

    cdheap.
      vector_push_back(SSHeapNode(cIndex,                    // cpointer is index
                     dIndexCard,                          // largest value of dsums <= target
                     c[cIndex].sum() + d[dIndexCard].sum())); // sum of heap element
  }

  cdheap.make_heap();
}

void initialize(const uint64_t S[MAXN], const int n,
                vector<SetNodeBitset> &a,vector<SetNodeBitset> &b,vector<SetNodeBitset> &c,vector<SetNodeBitset> &d,
                       SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper, const uint8_t maxCardinality) {


  initializeSets(S,n,a,b,c,d);
  initializeHeaps(a,b,c,d,abheap,cdheap,upper);
}


// ============================================================================
//
// Schroeppel and Shamir Helper Functions *** SUMS ***
//
// ============================================================================

void initializeSets(const uint64_t S[MAXN], const int n,
                    vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d) {

	generateAllSums (S, a, 0    , n/4-1);   // combinations of 1st quarter of numbers
	generateAllSums (S, b, n/4  , n/2-1); // combinations of 2nd quarter of numbers
	generateAllSums (S, c, n/2  , n*3/4-1); // combinations of 3rd quarter of numbers
	generateAllSums (S, d, n*3/4, n-1); // combinations of 4th quarter of numbers


}
void initializeHeaps(vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d,
                     SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper) {


	// --------------------------------------------------------------------------
	//construct min heap of initial combinations of sets from array A and array B
	// --------------------------------------------------------------------------

	for (size_t index = 0; index < a.size() && a[index] <= upper; index++) {  // for each set in A array
		abheap.vector_push_back(SSHeapNode(index,          // aptr pointer is same as index
													0,                    // smallest value of bsums is first element
													a[index]));           // a[index].sum + b.front().sum, b.front().sum = 0
		}
		abheap.make_heap();   // make the heap


		// --------------------------------------------------------------------------
		//construct max heap of initial combinations of sets from array C and array D
		// --------------------------------------------------------------------------

		int dIndex = d.size() - 1;                      // index to largest value in dsums

		for (size_t cIndex = 0; cIndex < c.size(); cIndex++) {  // process csums in increasing order
			if (c[cIndex] > upper) {
				break;                                    // no more entries to add to heap
			}

			while (c[cIndex] + d[dIndex] > upper) {
				dIndex--;
			}

			cdheap.
				vector_push_back(SSHeapNode(cIndex,                     // cpointer is index
											 dIndex,                  // largest value of dsums <= target
											 c[cIndex] + d[dIndex]));  // sum of heap element
		}
		cdheap.make_heap();
}

void initialize(const uint64_t S[MAXN], const int n,
                vector<uint64_t> &a,vector<uint64_t> &b,vector<uint64_t> &c,vector<uint64_t> &d,
                SSMinHeap &abheap, SSMaxHeap &cdheap, const uint64_t upper) {

  initializeSets(S,n,a,b,c,d);
  initializeHeaps(a,b,c,d,abheap,cdheap,upper);
}


// ============================================================================
//
// Schroeppel and Shamir Partition Algorithm
//
// ============================================================================

// Pop the top value from aHeap and replace with the proper next value
void nextAValue(SSMinHeap &aHeap, const vector<uint64_t> &a0, const vector<uint64_t> &a1) {

	if (aHeap.front().y() < a1.size() - 1) {                // there is another element in this column

		SSHeapNode newHeapNode(aHeap.front().x(),            // replace with next combination in this column
													 aHeap.front().y() + 1,
													 a0,a1);  // sum of new subset

		aHeap.pop();
		aHeap.push(newHeapNode);
	} else {                                                // Otherwise just pop node
		aHeap.pop();
	}
}

// Pop the top value from bHeap and replace with the proper next value
void nextBValue(SSMaxHeap &bHeap, const vector<uint64_t> &b0, const vector<uint64_t> &b1) {
	if (bHeap.front().y() > 0) {                         // there is another element in this column

	  			// Push new element
	  			bHeap.push(SSHeapNode(bHeap.front().x(),           // keep same row
										 	 	 	 	 	  bHeap.front().y() - 1,       // replace with next combination in column
										 	 	 	 	 	  b0,b1));                     // sum is c[x] + d[y]

	  			bHeap.pop();																			// Pop old one

	  		} else {                                              // no more elements in this column
	  			bHeap.pop();                                       // Remove old head
	  		}
}

uint64_t executeSS2(const partition::PartitionProblem &problem, ProblemStats &stats) {

	// *** FOR FORWARD ***
	//	const uint64_t *S = problem.S;

	// *** FOR REVERSE ***
	uint64_t *S = new uint64_t[problem.N];
	memcpy(S,problem.S,sizeof(uint64_t) * problem.N);
	// *** END FOR REVERSE ***


	const int N = problem.N;
	const uint64_t sum = problem.sum;
	uint64_t perfect = (sum+1) / 2;

  // Put first element in first subset
  uint64_t kkBound = kk(S, N, 2, sum);

  uint64_t best = getPartitionBest(kkBound, perfect, sum); 	// Get upper bound
  if (best == perfect) {  // If cgaBound was perfect, we are done
    return perfect;
  }
  uint64_t lower = getLowerBound(best,sum);											// Lower is complement of upper

  // Always put largest in set for efficiency. Since the first element has to
  // go in one set or the other, force it in
  uint64_t largest = S[0];


  // *** FOR REVERSE ***
  std::sort( S+1, S+N ,std::less<uint64_t>());  // Sort in ascending order
  // *** END FOR REVERSE ***

  vector<uint64_t> a0;    // all subset sums from first quarter of numbers
  vector<uint64_t> a1;    // all subset sums from second quarter of numbers
  vector<uint64_t> b0;    // all subset sums from third quarter of numbers
  vector<uint64_t> b1;    // all subset sums from fourth quarter of numbers

  SSMinHeap aHeap;  				// Min Heap Comparator
  SSMaxHeap bHeap;   				// Max Heap Comparator

  // Initialize the 4 vectors and 2 heaps, start from S+1 since
  // we always include the largest for efficiency
  initialize(S+1,N-1,a0,a1,b0,b1,aHeap,bHeap,best-1);

//  cout << endl << "a0: "; std::copy(a0.begin(),a0.end() , std::ostream_iterator<int>(std::cout," "));
//  cout << endl << "a1: "; std::copy(a1.begin(),a1.end() , std::ostream_iterator<int>(std::cout," "));
//  cout << endl << "b0: "; std::copy(b0.begin(),b0.end() , std::ostream_iterator<int>(std::cout," "));
//  cout << endl << "b1: "; std::copy(b1.begin(),b1.end() , std::ostream_iterator<int>(std::cout," "));
//  cout << endl << endl;
  uint64_t bValue = bHeap.front().sum();						// Current value from bHeap
  nextBValue(bHeap,b0,b1);	// Go to next b Value

  process_mem_usage(stats.residentMemory);

//  cout << "Perfect  : " << perfect << endl;
//  cout << "Initial B: " << bValue << endl;
  while (!aHeap.empty()) {                     // For each element of a
  	uint64_t aValue = aHeap.front().sum();		// Get top value from heap
  	uint64_t currentSum = largest + aValue + bValue;     //sum of elements on top of heaps
//  	if (currentSum < lower) {
//  	  cout << "A1: " << aValue << " B: " << bValue << " Sum: " << currentSum <<  " lb: " << lower << " ub: " << best << endl;
//  	}
  	while (currentSum >= lower && !bHeap.empty() ) {          // new set within bounds

//      cout << "A2: " << aValue << " B: " << bValue << " Sum: " << currentSum <<  " lb: " << lower << " ub: " << best << endl;

  		if (currentSum < best) {                               // if this value is better than best so far
  			best = getPartitionBest(currentSum, perfect, sum);    // Computer new upper bound from topsum
  			if (best == perfect) {                              	// If new topsum was perfect, we are done
  				return perfect;
  			}
  			lower = getLowerBound(best,sum);											// Lower is complement of upper

  		}

  		bValue = bHeap.front().sum();		// Get top value from heap
  		nextBValue(bHeap,b0,b1);				// Setup next B Value

  		currentSum = largest + aValue + bValue; // sum of two top elements of heaps
  	}

    if (currentSum >= lower) {	// This happens when we run out of elements in b, check remaining pairs with changing a
    	if (currentSum < best) {																	// If also <= upper, then we have a new bound
  			best = getPartitionBest(currentSum, perfect, sum);    	// Computer new upper bound from absum
  			if (best < perfect) {                              // If new absum was perfect, we are done
  				return perfect;
  			}
  			lower = getLowerBound(best,sum);											// Lower is complement of upper

    	} else {	// If abSum > upper, increasing a will not help
    		break;	// so we are done
    	}
    }


  	nextAValue(aHeap,a0,a1);
  }

  return best+1;   // Best is always 1 more that upper bound
}


} // end namespace
