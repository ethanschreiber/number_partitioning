/* This program implements the algorithm described in Mike Moffitt's
 IJCAI paper. This version runs KK to get the initial solution. In
 this version, SEARCH doesn't return the best solution cost, but
 simple sets the global variable BESTSOFAR to the cost of the best
 solution found so far. This version doesn't pass MAXSOFAR as an
 argument to the search routine, but maintains a global array of its
 value, indexed by k, the number of remaining subsets. */

#include "CGA_MW.hpp"

#include <stdio.h>                                    /* standard I/O library */
#include <math.h>                                     /* mathematical library */
#include <iostream>
#include <algorithm>
#include <vector>

namespace partition {

partition::CGA_MW::CGA_MW(const PartitionProblem& problem) :
		m_best(problem.sum), m_n(problem.N), m_k(problem.K), m_sum(problem.sum) {

	m_S = new uint64_t[problem.N];
	memcpy(m_S, problem.S, sizeof(uint64_t) * problem.N);
	std::sort(m_S, m_S + problem.N, std::greater<uint64_t>());
}

void CGA_MW::swapSort(uint64_t*& subsets, int idx) {
	for (int i=idx;i<m_k-1;i++) {
		if (subsets[i] > subsets[i+1]) {	// If out of order
			uint64_t tmp = subsets[i];		// swap i and i+1
			subsets[i] = subsets[i+1];
			subsets[i+1] = tmp;
		} else {							// keep going until in right spot
			break;
		}
	}
}

CGA_MW::~CGA_MW() {
	delete [] m_S;
}

uint64_t CGA_MW::search() {
	uint64_t *subsets = new uint64_t[m_k];
	memset(subsets,0,sizeof(uint64_t) * m_k);	// initialize to 0
	subsets[m_k-1] = m_S[0]; // place first integer in last subset

	uint64_t result = search(1, subsets, m_sum-m_S[0]);
	delete [] subsets;
	return result;
}


uint64_t partition::CGA_MW::search(int depth, uint64_t*& subsets, uint64_t sumRemaining) {

	uint64_t nextInteger = m_S[depth];
	uint64_t result = m_best;

//	cout << "Depth: " << depth << " Value: " << m_S[depth] << " Best: " << m_best << endl;

	if (depth == m_n-1) {					// if one integer left in S
//		for (int i=0;i<depth;i++) {
//					cout << "**";
//		}
//
//		cout << "Subsets: ";
//		for (int i=0;i<m_k;i++) {
//			cout << subsets[i] << " ";
//		}
//		cout << endl;

		uint64_t newValue = subsets[0] + nextInteger;
		uint64_t maxValue = std::max(subsets[m_k-1],newValue);
		m_best = std::min(m_best,maxValue); // Put value with smallest subset
		result = m_best;
	} else {								// If more than one left

		uint64_t newSumRemaining = sumRemaining - nextInteger;	// Remove new integer from sumRemaining remaining

		for (int i=0;i<m_k;i++) {                                   // For each subset

			if (i == 0 || subsets[i] != subsets[i-1]) {	// Don't consider both subsets if they have the same sum

				if (subsets[i] + nextInteger < m_best) {
					uint64_t *subsetsCopy = new uint64_t[m_k];			// Create copy of subsets
					memcpy(subsetsCopy,subsets,sizeof(uint64_t) * m_k);	// copy values
					subsetsCopy[i] += nextInteger;						// Add next value to proper subset

					swapSort(subsetsCopy,i);
//					for (int i=0;i<depth;i++) {
//						cout << "  ";
//					}

//					cout << "Added " << nextInteger << " to subset " << i << "  ";
//					cout << "Subsets: ";
//					for (int i=0;i<m_k;i++) {
//						cout << subsetsCopy[i] << " ";
//					}
//					cout << endl;

					result = std::min(m_best,search(depth+1,subsetsCopy,newSumRemaining));
					delete [] subsetsCopy;
				} else {
					break;
				}
			}
		}
	}
	return result;

}



uint64_t executeCGA_MW(const PartitionProblem &problem, ProblemStats &stats) {

	CGA_MW cga_mw(problem);
	return cga_mw.search();

} // end executeCGA_MW



} // end namespace


