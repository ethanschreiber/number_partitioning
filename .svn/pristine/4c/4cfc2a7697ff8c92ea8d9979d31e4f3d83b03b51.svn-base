/*
 * Horowitz_Sahni.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "Extended_Horowitz_Sahni.hpp"

#include <algorithm>

namespace ss {

size_t EHS(const uint64_t S[MAXN], const int n,
                      const uint64_t lower, const uint64_t upper,
                      vector<SetNodeBitset> &allSets)

{
  vector<SetNodeBitset> a;          // all subset sums from first half of numbers
  vector<SetNodeBitset> b;          // all subset sums from second half of numbers

  generateAllSets(S, n, a, 0    , n/2-1  ,true);      // combinations of 1st half of numbers (sorted ascending)
  generateAllSets(S, n, b, n/2  , n-1    ,false);     // combinations of 2nd quarter of numbers (sorted descending)

  size_t bAnchor = 0;                              // bPtr starts searching from here

  for (size_t aPtr=0; aPtr < a.size(); aPtr++) {
    const uint64_t &aValue = a[aPtr].sum();        // For readability, the value in the a vector

    while (bAnchor < b.size() &&                   // Move bAnchor to highest value that might fit
           aValue + b[bAnchor].sum() > upper) {    // Keep looking at smaller number in b until we
      bAnchor++;                                   // are under upper bound. (remember b is descending)
    }

    // Now add all aValue + b[bPtr] that are within range
    for (size_t bPtr=bAnchor;(bPtr < b.size()) && (aValue + b[bPtr].sum() >= lower); bPtr++) {
      allSets.push_back(SetNodeBitset(aValue + b[bPtr].sum(),a[aPtr].set() | b[bPtr].set()));
    }
  }

  return allSets.size();                                 // return the number of subset sums stored in ALLSUMS
}


size_t EHS(const uint64_t S[MAXN], const int n,
                      const uint64_t lower, const uint64_t upper,
                      vector<SetNodeVector> &allSets) {
  vector<SetNodeBitset> allSetsBitSet;

  EHS(S,n,lower,upper,allSetsBitSet);
  allSets.resize(allSetsBitSet.size());

  for (size_t i=0;i<allSetsBitSet.size();i++) {


    for (int j=0;j<n;j++) {
        if (allSetsBitSet[i].set()[j]) {
          allSets[i].push_back(S[j]);
        }
    }
  }

  return allSets.size();
}

}
