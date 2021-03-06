/*
 * BinPackingUtil.hpp
 *
 *  Created on: Nov 7, 2012
 *      Author: ethan
 */

#ifndef BINPACKINGUTIL_HPP_
#define BINPACKINGUTIL_HPP_

#include "PackingUtils.hpp"
#include <stdint.h>
namespace bp {
// ============================================================================
// TEST returns true if there is a subset of remaining included elements that
// when added to the subset sum so far is less than or equal to an excluded
// element, and for which the difference between the excluded element and the
// total subset sum is less than or equal to the residual capacity of the bin.
// Otherwise, it returns false.
// ============================================================================
bool testDominance(const uint64_t included[MAXN],   // array of included elements in descending order
          const int numInc,               // number of included elements
          const vector<uint64_t> &excluded,        // elements excluded that fit when they were excluded in descending order
          const int numExc,               // the number of such elements excluded
          const int excIdx,               // Initial index of current excluded element
          const uint64_t waste,            // Capacity left with all elements included
          const uint64_t sum);              // sum of elements already included in current subset

bool isDominated(
    const vector<uint64_t> &included, uint64_t includedSum,
    const vector<uint64_t> &excluded,
    const uint64_t capacity);    // Capacity left with first element included (rest are stored in included)

}// end namespace

// ============================================================================
// check to see if current set is a superset of any nogood. currentSet must end
// with -1, a sentinel value, before calling this.
// ============================================================================
bool checkNoGood(const SetNodeVector &currentSet, const vector<SetNodeVector> &noGoods,
                        const int binNum, vector<SetNodeVector> &newNoGoods,
                        const SetNodeVector bin[MAXN]);



bool processSolution(
		uint64_t &maxUsed,		  // The max capacity used in any bin (note, reference param)
		uint64_t &allowedWaste,     // max waste allowed. (note, reference param)
		int &bestSoFar,		      // The number of bins used in best solution so far (note, reference param)
		int N,					  // THe problem size
		uint64_t binCap,			  // The bin capacity
		int lowerBound,		      // lower bound on number of bins
		uint64_t sum,			  // The sum of all elements
		const SetNodeVector bin[MAXN],  // the set of elements in each bin
		const int binNum,         // index of the last bin filled + 1, starting from zero
		int solution[MAXN],      // bin index of each element in solution
		const uint64_t elements[MAXN],
		int minBins
		);



// ============================================================================
// REMOVEPAIRS - Removes all pairs which equal the capacity and stores the new
//               problem in processedProblem.
// ============================================================================

typedef vector<std::pair<uint64_t,uint64_t> > PairVector;
int removePairs(const BinPackingProblem &problem, BinPackingProblem &processedProblem,
                PairVector &perfectPairs);
#endif /* BINPACKINGUTIL_HPP_ */

