/*
 * HorowitzSahni.cpp
 *
 *  Created on: May 8, 2014
 *      Author: ethan
 */

#include "HorowitzSahni.hpp"
#include "SubsetSum.hpp"
#include "KK.hpp"
#include "../utils/MemoryUsage.hpp"
#include "../Utils.hpp"
namespace ss {

uint64_t executeHS2(const partition::PartitionProblem &problem, ProblemStats &stats) {

	const uint64_t *S = problem.S;
	const int N = problem.N;
	const uint64_t sum = problem.sum;

  // TODO: include largest element and search for perfect - largest
  // Will require some doing
  uint64_t perfect = (sum + 1) / 2;

  uint64_t kkBound = kk(S, N, 2, sum);

  uint64_t best = getPartitionBest(kkBound, perfect, sum); 	// Get best so far
  if (best == perfect) {  // If kkBound was perfect, we are done
    return perfect;
  }
  uint64_t lower = sum-(best-1);	// Lower is complement of upper

  // Always put largest in set for efficiency. Since the first element has to
  // go in one set or the other, force it in
  uint64_t largest = S[0];

  vector<uint64_t> a;                         // all subset sums from first half of numbers
  vector<uint64_t> b;                         // all subset sums from second half of numbers

  generateAllSums(S,  a, 1    , N/2   ,true);    // combinations of 1st half of numbers (sorted ascending) (skipping largest)
  generateAllSums(S,  b, N/2+1, N-1    ,false);  // combinations of 2nd half of numbers (sorted descending)

  auto bIt = b.begin();             // Iterator to b list
  uint64_t bValue = *bIt;						// Current value in blist
  bIt++;														// Set iterator to next value in b list

  process_mem_usage(stats.residentMemory);

  for (size_t i=0;i<a.size();i++) {
  	const uint64_t &aValue = a[i];
  	uint64_t currentSum = largest + aValue + bValue; // sum of a and b


    while (currentSum >= lower && bIt != b.end()) {								// While sum > lower bound, keep moving b pointer.

    	if (currentSum < best) {																	// If the abSum <= upper bound, compute new bounds
  			best = getPartitionBest(currentSum, perfect, sum);    	// Computer new upper bound from absum
  			if (best == perfect) {                              // If new absum was perfect, we are done
  				return perfect;
  			}
  			lower = sum-(best-1);                                  // Lower is complement of upper
    	}

    	bValue = *bIt;																				// Store the next b value
    	bIt++;																								// Increment the iterator
    	currentSum = largest + aValue + bValue;															// Compute the sum
    }

    if (currentSum >= lower) {	// This happens when we run out of elements in b, check remaining pairs with changing a
    	if (currentSum < best) {																	// If also <= upper, then we have a new bound
  			best = getPartitionBest(currentSum, perfect, sum);    	// Computer new upper bound from absum
  			if (best == perfect) {                              // If new absum was perfect, we are done
  				return perfect;
  			}
  			lower = sum-(best-1);                                // Lower is complement of upper

    	} else {	// If abSum > upper, increasing a will not help
    		break;	// so we are done
    	}
    }
  }
  return best+1;                                 // return the number of subset sums stored in ALLSUMS
}

}
