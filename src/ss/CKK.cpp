/*
 * SchroeppelShamirCompletion.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#include "CKK.hpp"

#include <algorithm>
#include <iostream>

#include <iomanip>
using std::cout;
using std::endl;
namespace ss {


// Merge two largest elements as difference, then copy the remaining
// from S into SCopy, inserting the difference in the correct spot
void diffFirstTwoAndCopy(const uint64_t S[], uint64_t SCopy[], const int n) {
  uint64_t diff2 = S[0] - S[1];  // difference of two largest elements

  // Merge two largest elements as difference, then copy the remaining
  // from into b, inserting the difference in the correct spot
  int i;                        // index into array
  for (i = 2; i < n; i++) {     // copy list and insert difference in order
    if (diff2 < S[i])
      SCopy[i - 2] = S[i];      // new number is less, keep copying
    else
      break;                    // found correct place, exit loop
  }
  SCopy[i - 2] = diff2;         // insert new element into array
  for (; i < n; i++) {     		// copy remaining elements
    SCopy[i - 1] = S[i];
  }
}

// Because we are using uint64_t, can't subtract first and
// then take absolute value
uint64_t absDiff(uint64_t a, uint64_t b) {
  if (a < b) {
    return b-a;
  } else {
    return a-b;
  }
}

/*
 * NOTE: S MUST be sorted in descending order when calling this.
 * CKK takes an array of uint64_t numbers A, its length N, and their TOTAL sum,
 and finds the best partition, leaving the resulting difference in the global
 variable ALPHA. It runs branch-and-bound, starting with the Karmarkar-Karp
 solution. At each point, the largest two remaining numbers are selected, and
 replaced with either their difference or their sum, representing assigning
 them to different sets or the same set, respectively. */

uint64_t ckk2(uint64_t S[],        // array of numbers
            const int n,            // number of elements in array
            const uint64_t sumRemaining,    // sum of all numbers
            uint64_t best)           // Don't set this when called externally
{
  // -----------------------------------
  // Base case when n ==4, KK is optimal
  // -----------------------------------
  if (n == 4) {
    uint64_t diffPart;           				         // Difference of partition
    uint64_t diff2 = S[0] - S[1];  			         // Difference of first two elements
    if (diff2 < S[2]) {							             // Compare diff to next largest
      diffPart = absDiff(S[2] - S[3], diff2);  // S[2] is biggest
    } else {
      diffPart = absDiff(diff2 - S[2], S[3]);  // diff2 is biggest
    }

    best = std::min(best, diffPart);		       // If better than best, then new best
  } else {

    // --------------------
    // Otherwise we recurse
    // --------------------

    { // *** Diff largest two ***
      uint64_t *SCopy = new uint64_t[n];  // new copy of list for recursive call
      diffFirstTwoAndCopy(S, SCopy, n); // Diff 2 largest and copy into SCopy in correct spot

      // We remove the first two and add back the diff of the first two:
      // newTotal = total - sum[0] - sum[1] + (sum[0]-sum[1])
      //          = total - sum[0] + sum[0] - sum[1] - sum[1]
      //          = total - sum[1] - sum[1]
      uint64_t newSumRemaining = sumRemaining - S[1] - S[1];

      // compare the first element to the sum of the rest:
      const uint64_t &largest = SCopy[0];
      uint64_t        sumRest = newSumRemaining - SCopy[0]; // sum of all elements except first
      if (largest >= sumRest) {             // if largest element >= sum of rest
        delete [] SCopy;
        return std::min(best, largest - sumRest);

      } else {    // if difference is larger than rest, so is the sum

        best = ckk2(SCopy, n - 1, newSumRemaining,best);            // one less element, new total
        delete [] SCopy;
      }
    }

    // If the diff recursion didn't find a solution
    if (best > 1)
    { // compare the first element to the sum of the rest:
      uint64_t largest = S[0] + S[1];          // sum of largest two numbers is new largest
      uint64_t sumButTwo = sumRemaining - S[0] - S[1];  // sum of all elements except first two
      if (largest >= sumButTwo) {               // if largest element >= sum of rest
        return std::min(best, largest - sumButTwo);
      }

      S[1] = largest;              // sum of two largest is new largest element
      best = ckk2(S + 1, n - 1, sumRemaining, best);    // call on subarray with one less element
      S[1] = largest - S[0];       // restore array to previous state
    }
  }
  return best;
}


uint64_t executeCKK2(uint64_t S[],        // array of numbers
                    const int n,            // number of elements in array
                    const uint64_t sum) {    // sum of all numbers

  uint64_t diff = ckk2(S,n,sum);

  return (sum + diff) / 2;      // Return larger of two sets
}
} /* namespace ss */

//const int SIZE = 40;
//const int NUM_TRIALS = 100;
//int main() {
//  uint64_t S[SIZE];
//
//  for (int trial = 0; trial < NUM_TRIALS; trial++) {
//
//    uint64_t sum = 0;
//    for (int i = 0; i < SIZE; i++) {
//      S[i] = rand();
//      sum += S[i];
//    }
//    std::sort(S,S+SIZE,std::greater<uint64_t>());
//    cout << "Trial " << std::setw(3) << trial << ": ";
//
//    uint64_t alpha1 = ss::ckk(S, SIZE, sum,1);
//
//    cout << "Alpha1: " << alpha1 << endl;
//
//  }
//
//  return 0;
//}
