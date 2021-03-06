/*
 * DynamicProgramming.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "DPPartition.hpp"
#include "../utils/MemoryUsage.hpp"
#include <iostream>
#include <cstring> // for memcpy
#include <iomanip>
using std::vector;
using std::cout;
using std::endl;
namespace ss {


uint64_t executeDPPartition(const partition::PartitionProblem &problem, ProblemStats &stats) {

	const uint64_t *SInit = problem.S;
	const int n = problem.N;
	const uint64_t sum = problem.sum;

  uint64_t *S = new uint64_t[n];              	// A copy of SInit
  memcpy(S,SInit,sizeof(uint64_t) * n);					// Copy the values
  std::sort(S,S+n, std::greater<uint64_t>()); 	// Sort m_S descending

  uint64_t best = sum;                        	// The best value so far
  const uint64_t PERFECT = (sum+1) / 2;       	// The value of a perfect partition
  const uint64_t MIN_INPUT_VALUE = S[n-1];      			// Last element is smallest
  const size_t ROW_SIZE = PERFECT - MIN_INPUT_VALUE;	// Need one bit for each value between m_smallestValue and PERFECT

  // Create bitset for the dynamic programming. Represents all subsets in descending order
  // from PERFECT-1 at idx 0 to MIN_VALUE at idx ROW_SIZE-1.
  DPBitset row(ROW_SIZE);

  uint64_t MIN_INDEX = ROW_SIZE;
  uint64_t MAX_INDEX = 0;

  for (int i=0;i<n;i++) {                         									// For each integer in S
    uint64_t &s_i = S[i];																						// Get reference to integer

		size_t idx = MIN_INDEX;                													// Iterate through all previously set values
		while (idx <= MAX_INDEX && idx != DPBitset::npos) {             // While there are more values left

			uint64_t value = idxValue(ROW_SIZE, MIN_INPUT_VALUE, idx);		// get value corresponding to idx
			uint64_t newValue = value + s_i;															// Add new value

			if (newValue < PERFECT) {            									       	// If less than PERFECT

				uint64_t newIdx = idxValue(ROW_SIZE, MIN_INPUT_VALUE, newValue);	// Get new Idx corresponding to new value
				row.set(newIdx,true); 																			// Set new bit in row at newIdx
				MIN_INDEX = std::min(MIN_INDEX, newIdx);
				MAX_INDEX = std::max(MAX_INDEX, newIdx);
			  newValue = sum - newValue;                									// use complement to compare to best
			} else if (newValue == PERFECT) {            									// Check for perfect partition
			  best = newValue;
				goto CLEANUP;																								// If perfect, cleanup and return
			} else {
				best = std::min(newValue,best);           									// keep track of best so far
			}

		  idx = row.find_next(idx);																			// Go to next idx
		}

		// Put in the new lone element
		idx = idxValue(ROW_SIZE,MIN_INPUT_VALUE,s_i);										// idx of s_i
		row.set(idx,true);					  																	// Set bit for s_i
		MIN_INDEX = std::min(MIN_INDEX, idx);
		MAX_INDEX = std::max(MAX_INDEX, idx);

	}

  CLEANUP:
  process_mem_usage(stats.residentMemory);
  delete [] S;
	return best;
}

} // end namespace
