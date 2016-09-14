// This program takes a random set of N 48-bit integers, and lower and
//   upper bounds on subset sums, and generates all possible sums within
//   the bounds of subsets of the given numbers.  It is an extension of
//   the Schroeppel and Shamir algorithm.

#include "SubsetSum.hpp"
#include <stdio.h>                                    // standard I/O library
#include <stdint.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>


#include <iostream>

using std::cout;

using std::endl;

namespace ss {

// This is the max cardinality version.
void generateAllSetsRecursive (const uint64_t S[],
                      vector<SetNodeBitset> &sets,
                      int first, int last,
                      uint64_t curSum, DynamicBitset &curSet,
                      const uint8_t card, const uint8_t maxCard )

{
  if (first > last) {                                    // set is completed
    sets.push_back(SetNodeBitset(curSum,curSet));
  } else {                                               // set not yet completed

    // exclude next element
    generateAllSetsRecursive (S, sets, first+1, last, curSum, curSet, card, maxCard);

    if (card < maxCard) {
      // Include next element
      curSet[first] = true;
      generateAllSetsRecursive (S, sets, first+1, last, curSum + S[first], curSet, card+1, maxCard);
      curSet[first] = false;
    }
  }
}


// Helper function, provides missing variables
// also sorts sets before returning
void generateAllSets(const uint64_t S[], const int N,
             	 	 	 	 vector<SetNodeBitset> &sets,
             	 	 	 	 int first, int last,
             	 	 	 	 bool sortAscending, uint8_t maxCardinality) {
  DynamicBitset curSet(N);                                   // The number of total input elements
  generateAllSetsRecursive(S,sets,first,last,0ll,curSet,0,maxCardinality);
  if (sortAscending) {
    std::sort( sets.begin(), sets.end() ,std::less<SetNodeBitset>());  // Sort in ascending order
  } else {
    std::sort( sets.begin(), sets.end() ,std::greater<SetNodeBitset>());  // Sort in descending order
  }

}



// GENSUMS takes an array NUMS of numbers, the index of the LAST
//     number, and generates an array SUMS of all sums that can be
//     achieved by adding together numbers from NUMS.  Each number in NUMS
//     can only be used at most once.  To allow the recursion, additional
//     arguments include an index NEXT to the next element of the NUMS
//     array, and the SUMSOFAR achieved down this path.  There is also a
//     global pointer NEXTSUM to the next empty element of the SUMS
//     array. At the end, it is equal to the number of sums created.

void generateAllSumsRecursive (const uint64_t nums[],         // array of original numbers
              vector<uint64_t> &sums,                // array of sums
              int next,                              // pointer to next element of NUMS array
              int last,															 // index of last element in array
              uint64_t curSum) {                     // sum so far of elements in current subset


  if (next == last) {                                // reached last element of array
    sums.push_back(curSum);                        // don't add last element
    sums.push_back(curSum + nums[next]);           // add last element to sum
  } else {
    generateAllSumsRecursive (nums, sums, next+1, last, curSum);   // don't add last element
    generateAllSumsRecursive (nums, sums, next+1, last, curSum + nums[next]); // add last element
  }
}

// Helper function, provides missing variables
// also sorts sets before returning
void generateAllSums(const uint64_t S[],
             	 	 	 	 vector<uint64_t> &sums,
             	 	 	 	 int first, int last,
             	 	 	 	 bool sortAscending) {
  generateAllSumsRecursive(S,sums,first,last,0ll);
  if (sortAscending) {
    std::sort( sums.begin(), sums.end() ,std::less<uint64_t>());  // Sort in ascending order
  } else {
    std::sort( sums.begin(), sums.end() ,std::greater<uint64_t>());  // Sort in descending order
  }
}



} // end namespace
