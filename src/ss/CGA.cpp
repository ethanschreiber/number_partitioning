/* This program finds optimal solutions to arbitrary multi-way number
 partitioning problems using the complete greedy algorithm.  It sorts
 the numbers, and assigns each to one of the sets in turn, searching
 the resulting k-ary tree, using branch-and-bound.  The first solution
 found is the greedy solution.  This version checks to see if any of
 the current subsets are empty, and if so, only places the next
 element in one empty subset. This version uses all uint64_t
 variables, except for indices.  This version uses as the optimization
 criteria the maximum subset sum instead of the difference between the
 largest and smallest subset sums.  This version is optimized to
 eliminate all but one permutation of the subsets. Rather than making
 a new copy of the array of subset sums for each recursive call, this
 version keeps a single global array, makes changes at each node, and
 undoes those changes for the next node. This version prevents
 assigning a number to multiple empty sets by passing as a parameter
 the first subset to assign a new number to.  This version computes
 the greedy solution when there are two numbers left, instead of
 waiting until there is only one number left. This version uses 48-bit
 integers. */

#include "CGA.hpp"
#include <stdio.h>            // standard I/O library
#include <stdint.h>
#include <cstring>            // for memset
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <queue>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

namespace ss {

static uint64_t g_nodeCount;          // total calls to search routine

// ==========================================================
// Complete Greedy Algorithm for the 2-way partition problem.
// ==========================================================
uint64_t cga2 (const uint64_t *S, const size_t N, const uint64_t sumRemaining,
								 const uint64_t smallSum, const uint64_t largeSum,
								 const uint64_t perfect, const uint64_t ub) {

	g_nodeCount++;

	uint64_t returnValue=0;
	if (N == 1) {																// If one element left
		returnValue = std::max(smallSum+S[0],largeSum); 								// Put it in smaller set

	} else  if (sumRemaining + smallSum < largeSum) {	// If all remaining fit in smaller sum, prune
		returnValue = largeSum;
	} else {

	  // Put the new integer S[0] in the smaller subset
	  returnValue = cga2(S+1, N-1, sumRemaining - S[0],
	  									 std::min(smallSum+S[0], largeSum),    // The smaller of smallSubset+S[0] and largeSubset
	  									 std::max(smallSum+S[0], largeSum),    // The larger of smallSubset+S[0] and largeSubset
	  									 perfect,ub);

		// If the subset sums are equal, only need to add element to the first one
		if (smallSum != largeSum && returnValue > perfect) {

		  // Put the new integer S[0] in the larger subset
		  // largeSubset + S[0] is always greater than smallSubset
			returnValue = std::min(returnValue,
														cga2(S+1, N-1, sumRemaining - S[0], smallSum, largeSum+S[0], perfect,ub));
		}
	}
	return returnValue;

}

uint64_t executeCGA2(const uint64_t S[], const int N, const uint64_t sum) {

	uint64_t perfect = (sum+1) / 2;
	g_nodeCount = 0;
	// Put first element in first subset
	uint64_t solution = cga2(S+1,N-1, sum-S[0], 0, S[0], perfect, sum);

	//cout << "Node Count: " << g_nodeCount << endl;
	return solution;
}


// =============================================================
// The greedy heuristic, assumes S is sorted in descending order
// =============================================================

uint64_t printQueue(std::priority_queue<uint64_t,vector<int>,std::greater<int> > q) {
  while (!q.empty()) {
    cout << q.top() << " ";
    q.pop();
  }
  cout << endl;
}

uint64_t executeGreedy(const uint64_t S[], const int N, const int K) {

  std::priority_queue<uint64_t,vector<int>,std::greater<int> > q;		// Min Heap
  uint64_t maxValue = S[0];																					// Max starts as largest element
  for (int i=0;i<K;i++) {		// Put first K on heap
  	q.push(S[i]);

//  	cout << "S[" << i << "]: " << std::setw(2) << S[i] << " - ";printQueue(q);
  }


  // Iterate through rest
  for (int i=K;i<N;i++) {
  	uint64_t newValue = q.top() + S[i];				// New value is min value plus next integer of S
  	maxValue = std::max(maxValue,newValue);		// If greater than maxValue, update it
  	q.pop();																	// Remove old value
  	q.push(newValue);													// Push back new value
//  	cout << "S[" << i << "]: " << std::setw(2) << S[i] << " - "; printQueue(q);
  }

  return maxValue;
}

} // end namespace
