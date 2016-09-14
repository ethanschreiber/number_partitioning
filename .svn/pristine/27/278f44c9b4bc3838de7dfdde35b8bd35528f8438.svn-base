/*
 * BinPackingUtil.cpp
 *
 *  Created on: Nov 7, 2012
 *      Author: ethan
 */

#include "BinCompletionUtils.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include "completions/IECompletionGenerator.hpp"
using std::vector;
using std::cout;
using std::endl;

// ============================================================================
// testDominance returns true if there is a subset of remaining included
// elements that when added to the subset sum so far is less than or equal to
// an excluded element, and for which the difference between the excluded
// element and the total subset sum is less than or equal to the wasted
// capacity of the bin. Otherwise, it returns false.
// ============================================================================


//if (testDominance(&included[firstIdx], numInc - firstIdx, excluded, numExc,
//                  numExc - 1, residual - included.sum, 0)) {

  bool bp::testDominance(const uint64_t included[MAXN],   // array of included elements in descending order
          const int numInc,               // number of included elements
          const vector<uint64_t> &excluded,        // elements excluded that fit when they were excluded in descending order
          const int numExc,               // the number of such elements excluded
          const int excIdx,               // Initial index of current excluded element
          const uint64_t waste,           // Capacity left with all elements included
          const uint64_t sum)             // sum of elements already included in current subset
{
  if (numInc >= 2) {                          // there are at least two numbers remaining
    if (testDominance(&included[1], numInc-1, excluded, numExc, excIdx, waste, sum)) {
      return true;                            // exclude first element
    }
  }
  uint64_t newSum = sum + included[0];         // add first included element

  //cout << "New Sum: " << newSum << endl;
  if (newSum > excluded[0]) { // If sum larger than all elements
    return false;              // Then no dominance possible
  }

  int newExcIdx;
  for (newExcIdx = excIdx; newExcIdx >= 0; newExcIdx--) { // consider excluded elements in increasing order of size
    if (newSum <= excluded[newExcIdx] && included[0] != excluded[newExcIdx]) {
      break;                                              // first excluded element >= sum
    }
  }


  if (newExcIdx < 0) {
    return false;                                // no such excluded element, no dominance
  }


  if (excluded[newExcIdx] - newSum <= waste &&   //  subset with excluded element won't exceed capacity
      excluded[newExcIdx] != included[0]) {      // and excluded element is not a copy of included element
    return true;                                 // then this set is dominated
  }

  if (numInc > 1) {                             // there are more elements to include
    if (testDominance(&included[1], numInc-1, excluded, numExc, newExcIdx, waste, newSum)) {
      return true;                              // include next element
    }
  }

  return false;                                 // didn't find any dominance, return no dominance
}


bool bp::isDominated(
    const vector<uint64_t> &included, uint64_t includedSum,
    const vector<uint64_t> &excluded,
    const uint64_t waste)    // Capacity left with first element included (rest are stored in included)
{

  int numExc = excluded.size();
  if (numExc > 0) { // if any excluded elements, test for dominance

    int numInc = included.size();
    int firstIdx = 0; // index of first included element <= largest excluded element

    // find index of first included element <= biggest excluded element
    while (firstIdx < numInc && included[firstIdx] > excluded[0]) {
      firstIdx++;
    }

    if (firstIdx < numInc) {                       // there is an included element smaller than biggest excluded

      if (testDominance(&included[firstIdx], numInc - firstIdx, excluded, numExc,
                        numExc - 1, waste, 0)) {
        return true;
      }

    }
  }

  return false;
}


// ============================================================================
// check to see if current set is a superset of any nogood. currentSet must end
// with -1, a sentinel value, before calling this.
// ============================================================================
bool checkNoGood(const SetNodeVector &currentSet, const vector<SetNodeVector> &noGoods,
                        const int binNum, vector<SetNodeVector> &newNoGoods,
                        const SetNodeVector bin[MAXN]) {

  for (size_t i = 0; i < noGoods.size(); i++) { // for each existing nogood

    int index = 0;                              // index of current member of thiset
    bool isSubset = true;                       // initially, nogood is a subset of completion

    //matches first element of bin,
    // true if thiset and nogood have an element in common, false otherwise
    bool common = (noGoods[i].member[0] == bin[binNum].member[0]);

    int nogoodCard = noGoods[i].cardinality();
    for (int j = (common ? 1 : 0); j < nogoodCard; j++) {           //for each remaining element of nogood set
      while (currentSet.member[index] > noGoods[i].member[j]) {
        index++; // look for element
      }
      if (currentSet.member[index] < noGoods[i].member[j]) {  //element not in completion
        isSubset = false;                                     // nogood is not a subset of completion
        if (common) {
          goto OVERLAP;
        }
      } else {                              // found matching element
        common = true;                      // sets have at least one element in common
        if (!isSubset) {
          goto OVERLAP;                     // but nogood is not a subset
        }
        index++;
      }
    }                                       // advance to next element of completion

    if (isSubset) {
      return true;                          // this nogood is a subset of this completion
    }

    newNoGoods.push_back(noGoods[i]);       // nogood and completion have no overlap

    OVERLAP: continue;
  } // nogood and thiset overlap, but nogood not a subset

  return false;
}

// ============================================================================
// This processes the solution.  If it is the best so far, it translates the
// sets of elements into the solution vector, where each element contains the
// index of the bin that contains it. It also updates the allowedwaste in that case.
// It returns true if the solution is proven optimal, false otherwise
// ============================================================================
bool processSolution(
		uint64_t &maxUsed,		         // The max capacity used in any bin (note, reference param)
		uint64_t &allowedWaste,        // max waste allowed. (note, reference param)
		int &bestSoFar,		             // The number of bins used in best solution so far (note, reference param)
		int N,					               // THe problem size
		uint64_t binCap,			         // The bin capacity
		int lowerBound,		             // lower bound on number of bins
		uint64_t sum,			             // The sum of all elements
		const SetNodeVector bin[MAXN],       // the set of elements in each bin
		const int binNum,              // index of the last bin filled + 1, starting from zero
		int solution[MAXN],            // bin index of each element in solution
		const uint64_t elements[MAXN], //
		int minBins)
{
  uint64_t elementsCopy[MAXN];     // copy of original elements

  // either (a) ran out of elements before bins or (b) processed at least all initial bins
  int numBins = (minBins > binNum) ? minBins : binNum;

  //cout << "   Num Bins: " << numBins << endl;
  assert(numBins >0 && "Need at least one bin!");

  if (numBins < bestSoFar) {                             // solution uses fewer bins than best so far
    memcpy(elementsCopy,elements,sizeof(uint64_t) * N);  // Copy elements

    for (int i = 0; i < numBins; i++) {                  // for each bin index

      int cardinality = bin[i].cardinality();
      uint64_t sum = 0;
      for (int j = 0; j < cardinality; j++) {            // for each element of the set

        uint64_t target = bin[i].member[j];              // the numeric element in the bin
        sum += target;
        int k;
        for (k = 0; k < N; k++) {                        // look for index of element in array
          if (elementsCopy[k] == target) {
            break;
          }
        }
        solution[k] = i;                                // assign the bin index to the solution array
        elementsCopy[k] = -1;                           // scratch this copy out as used
      }
      maxUsed = std::max(maxUsed,sum);                  // Set max used to largest capacity used in any bin
    }
    bestSoFar = numBins;                                // update bestsofar to new best value
    allowedWaste = (bestSoFar-1) * binCap - sum;        // maximum waste to do better

  }

  // We use <= instead of == because lowerBound can be specified as a "good enough" value
  // instead of a true lowerbound. This is useful if using this for number partitioning. Ir you
  // want an optimal packing, make sure lowerBound is a true lowerBound for the problem.
  return (bestSoFar <= lowerBound);                     // true if optimal, false otherwise
}


// ============================================================================
// REMOVEPAIRS - Removes all pairs which equal the capacity and stores the new
//               problem in processedProblem. Adds these pairs to perfectPairs.
// ============================================================================

int removePairs(const BinPackingProblem &problem, BinPackingProblem &processedProblem,
                PairVector &perfectPairs) {

  int pairsRemoved = 0;
  int N = problem.N;
  uint64_t S[MAXN];

  int headIdx = 0;
  int tailIdx = problem.N - 1;

  // Copy numbers from problem.S to S. If two numbers sum to capacity, set the numbers
  // to -1 to be removed next
  while (headIdx < tailIdx) {
    uint64_t pairSum = problem.S[headIdx] + problem.S[tailIdx];
    if (pairSum < problem.capacity) { // If too small, move the tail idx
      S[tailIdx] = problem.S[tailIdx];
      tailIdx--;
    } else if (pairSum > problem.capacity) {// If too large, move the head idx
      S[headIdx] = problem.S[headIdx];
      headIdx++;
    } else {          // If just right, remove both head and tail
      perfectPairs.push_back(std::pair<uint64_t,uint64_t>(problem.S[headIdx],problem.S[tailIdx]));
      S[headIdx] = -1;
      S[tailIdx] = -1;
      headIdx++;
      tailIdx--;
      N -= 2;
      pairsRemoved++;
    }
  }

  // Copy last number
  if (headIdx == tailIdx) {
    S[headIdx] = problem.S[headIdx];
  }

  // Remove -1's
  int toIdx = 0;
  for (int i = 0; i < problem.N; i++) {
    if (S[i] > 0) {
      S[toIdx] = S[i];
      toIdx++;
    }
  }

  processedProblem.reset(S,N);

  return pairsRemoved;
}
