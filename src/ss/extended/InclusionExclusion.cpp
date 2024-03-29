/*
 * InclusionExclusion.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */


#include "InclusionExclusion.hpp"

#include <algorithm>


namespace ss {


//size_t generateSetsIE(const uint64_t S[MAXN], const int SIdx, const int n,
//                      const uint64_t sum, const uint64_t sumRemaining, DynamicBitset &bitset,
//                      const uint64_t lower, const uint64_t upper,
//                      vector<SSNode> &allSets) {
//  size_t count = 0;
//
//  // If it is possible to find more subset sums in range, recurse
//  if (sum + sumRemaining >= lower && (SIdx < n)) {
//
//    // Exclude first
//    count += generateSetsIE(S, SIdx+1, n, sum          , sumRemaining - S[SIdx],
//                            bitset, lower, upper, allSets);
//
//    // Now include
//    uint64_t newSum = sum + S[SIdx];   // The new sum with the included element
//
//    if (newSum <= upper) {
//      bitset[SIdx] = true;            // Set corresponding bit in bitset
//
//      if (newSum >= lower) {          // In range, increase count
//        allSets.push_back(SSNode(sum,bitset));
//        count++;
//      }
//
//      // Recurse
//      count += generateSetsIE(S, SIdx+1, n, sum + S[SIdx], sumRemaining - S[SIdx],
//                              bitset, lower, upper, allSets);
//      bitset[SIdx] = false;           // Unset bit
//    }
//  }
//
//  return count;
//}


size_t incExc(const uint64_t S[MAXN], const int SIdx, const int n,
                      const uint64_t sum, const uint64_t sumRemaining, DynamicBitset &bitset,
                      const uint64_t lower, const uint64_t upper,
                      vector<SetNodeBitset> &allSets) {
  size_t count = 0;

  // If it is possible to find more subset sums in range, recurse
  if (sum + sumRemaining >= lower && (SIdx < n)) {

    // Now include
    uint64_t newSum = sum + S[SIdx];   // The new sum with the included element

    if (newSum <= upper) {


      bitset[SIdx] = true;            // Set corresponding bit in bitset

      if (newSum >= lower) {          // In range, increase count
        allSets.push_back(SetNodeBitset(sum,bitset));
        count++;
      }

      // Recurse
      count += incExc(S, SIdx+1, n, sum + S[SIdx], sumRemaining - S[SIdx],
                              bitset, lower, upper, allSets);
      bitset[SIdx] = false;           // Unset bit


    }
      // Exclude
      count += incExc(S, SIdx+1, n, sum          , sumRemaining - S[SIdx],
                              bitset, lower, upper, allSets);


  }

  return count;
}



size_t incExc(const uint64_t S[MAXN], const int n, const uint64_t sumRemaining,
                      const uint64_t lower, const uint64_t upper,
                      vector<SetNodeBitset> &allSets)
{

  DynamicBitset bitset(n);
  if (lower == 0) {
    allSets.push_back(SetNodeBitset(0,bitset));
    return 1 + incExc(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);
  } else {
    size_t x = incExc(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);;

    return  x;
  }
}



// ================
// Simple dominance
// ================
size_t generateSetsIESimpleDominance(const uint64_t S[MAXN], const int SIdx, const int n,
                                     const uint64_t sum, const uint64_t sumRemaining, DynamicBitset &bitset,
                                     const uint64_t lower, const uint64_t upper,
                                     vector<SetNodeBitset> &allSets) {
  size_t count = 0;

  if (SIdx == n) {                                    // If at a leaf
    if (sum >= lower && sum <= upper) {               // If in range
        allSets.push_back(SetNodeBitset(sum,bitset));     // Add set
        count++;                                      // 1 more added
    }
  } else if (sum + sumRemaining >= lower && (SIdx < n)) {   // If still possible to get in range

    if (sum + sumRemaining <= upper) {
      for (int i=SIdx;i<n;i++) {
        bitset[i] = true;
      }
      allSets.push_back(SetNodeBitset(sum+sumRemaining,bitset));     // Add set
      count++;                                      // 1 more added

      for (int i=SIdx;i<n;i++) {
        bitset[i] = false;
      }
    } else {

      // First Include
      if (sum + S[SIdx] <= upper) {     //

        bitset[SIdx] = true;            // Set corresponding bit in bitset

        count += generateSetsIESimpleDominance(S, SIdx+1, n, sum + S[SIdx], sumRemaining - S[SIdx],
                                               bitset, lower, upper, allSets);
        bitset[SIdx] = false;           // Unset bit
      }

        count += generateSetsIESimpleDominance(S, SIdx+1, n, sum          , sumRemaining - S[SIdx],
                                               bitset, lower, upper, allSets);
    }

  }

  return count;
}



size_t generateSetsIESimpleDominance(const uint64_t S[MAXN], const int n, const uint64_t sumRemaining,
                                     const uint64_t lower, const uint64_t upper,
                                     vector<SetNodeBitset> &allSets)
{

  DynamicBitset bitset(n);
  if (lower == 0) {
    allSets.push_back(SetNodeBitset(0,bitset));
    return 1 + generateSetsIESimpleDominance(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);
  } else {
    size_t x = generateSetsIESimpleDominance(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);;

    return  x;
  }
}


// =================
// Moffitt dominance
// =================
size_t generateSetsIEMoffittDominance(const uint64_t S[MAXN], const int SIdx, const int n,
                                     const uint64_t sum, const uint64_t sumRemaining, DynamicBitset &bitset,
                                     const uint64_t lower, const uint64_t upper,
                                     vector<SetNodeBitset> &allSets) {
  size_t count = 0;

  if (SIdx == n) {                                    // If at a leaf
    if (sum >= lower && sum <= upper) {               // If in range
        allSets.push_back(SetNodeBitset(sum,bitset));     // Add set
        count++;                                      // 1 more added
    }
  } else if (sum + sumRemaining >= lower && (SIdx < n)) {   // If still possible to get in range

    if (sum + sumRemaining <= upper) {
      for (int i=SIdx;i<n;i++) {
        bitset[i] = true;
      }
      allSets.push_back(SetNodeBitset(sum+sumRemaining,bitset));     // Add set
      count++;                                      // 1 more added

      for (int i=SIdx;i<n;i++) {
        bitset[i] = false;
      }
    } else {

      // First Include
      if (sum + S[SIdx] <= upper) {     //

        bitset[SIdx] = true;            // Set corresponding bit in bitset

        count += generateSetsIEMoffittDominance(S, SIdx+1, n, sum + S[SIdx], sumRemaining - S[SIdx],
                                               bitset, lower, upper, allSets);
        bitset[SIdx] = false;           // Unset bit

        // Now exclude
        count += generateSetsIEMoffittDominance(S, SIdx+1, n, sum          , sumRemaining - S[SIdx],
                                                bitset,
                                                std::max(lower,sum + S[SIdx] + 1 ), upper, allSets);
      } else {
        // Now exclude
        count += generateSetsIEMoffittDominance(S, SIdx+1, n, sum          , sumRemaining - S[SIdx],
                                                bitset, lower, upper, allSets);
      }


    }

  }

  return count;
}



size_t generateSetsIEMoffittDominance(const uint64_t S[MAXN], const int n, const uint64_t sumRemaining,
                                     const uint64_t lower, const uint64_t upper,
                                     vector<SetNodeBitset> &allSets)
{

  DynamicBitset bitset(n);
  if (lower == 0) {
    allSets.push_back(SetNodeBitset(0,bitset));
    return 1 + generateSetsIEMoffittDominance(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);
  } else {
    size_t x = generateSetsIEMoffittDominance(S,0,n,0,sumRemaining,bitset,lower,upper,allSets);;

    return  x;
  }
}


} // end namespace
