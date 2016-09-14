/*
 * OrderedPowerSet.cpp

 Here's an algorithm. The basic idea is that each number in the original set
 iterates through the list of subsets you've already found, trying to see if
 adding that number to the subset it's currently considering results in the
 smallest subset sum not yet found.

 The algorithm uses four arrays (all of which are indexed starting with 0).

 S consists of the numbers in the original set; i.e., S=[1,4,5,9] in your example.
 L is the list of subsets found so far.
 A[i] contains the subset that S[i] is currently considering.
 SUM[i] is the sum of the elements of subset i in L.

 Algorithm:

 1. Initialize S to numbers in the original set, all entries of A to 0, L[0]={},
 SUM[0]=0. Let j=1.
 2. For iteration j find the minimum of SUM[A[i]]+S[i] over all numbers S[i] in
 the original set. (This finds the subset with smallest sum not yet in L.)
 Tie-breaking is done by number of elements in the set. Let i∗ denote the
 argmin.

 3. Let L[j]=L[A[i∗]]∪{S[i∗]}. Let SUM[j]=SUM[A[i∗]]+S[i∗]. (This updates L and SUM
 with the new subset.)

 4. Increase A[i∗] to the next item in L that has no number larger than S[i∗].
 If there is none, let A[i∗]= NULL. (This finds the next subset in L to
 consider for the number S[i∗] just added to an existing subset in L to
 create the subset just added to L.)

 5. If all entries in A[i] are NULL, then stop, else increment j and go to
 Step 2.

 *
 *  Created on: Jun 11, 2013
 *      Author: ethan
 */

#include "OrderedPowerSet.hpp"
#include "SubsetSum.hpp"
#include <iostream>
#include <deque>
#include <iomanip>
#include <cstring> // for memset
using std::cout;
using std::endl;
using std::deque;


namespace ss {

size_t generateSetsOPS(const uint64_t S[MAXN], const int N,
    const uint64_t lower, const uint64_t upper, deque<SetNodeBitset> &LParam) {

  PersisistentIndexDeque<SetNodeBitset> L(LParam);   // Deque whose indexing stays the same even with pop_front() is called
  L.push_back(SetNodeBitset(N));                     // Add empty subset

  OPSMinHeap nextHeap;                        // Next element to add

  // Initialize the heap
  for (int i=0;i<N;i++) {                     // For each input element
    nextHeap.push(OPSHeapNode(0,i,S[i],0));   // One entry in heap pointing to 0 idx of L, i idx of S, sum of S[i] and cardinality 0
  }

  do {
    OPSHeapNode node = nextHeap.front();      // Get min sum
    nextHeap.pop();                           // Pop it from heap

    if (node.sum <= upper &&                  // Add if less than upper bound AND
        (node.sum >= lower ||                 // greater than upper bound OR
         node.SIdx < (size_t) N-1)) {         // another set can be created from this one.
      L.push_back(SetNodeBitset(L[node.LIdx]));      // 3. Insert copy of old entry into L
      L.back().set(node.SIdx, S[node.SIdx]);  // Insert new element into copy (changes sum too)
    }

    if (node.SIdx == (size_t) N-1 &&     // If this node corresponds to the largest value in S
        L[node.LIdx].sum() < lower) {    // And the value is less than the lower bound

      L.pop_front();
    }

    size_t LIdx = node.LIdx;
    do {                                  // Find next set we can add S[node.SIdx] to.
      LIdx++;
    } while (LIdx < L.pushCount() && L[LIdx].max() >= node.SIdx);

    if (LIdx < L.pushCount()) {      // If there is another subset to point to
      nextHeap.push(OPSHeapNode(LIdx,node.SIdx,L[LIdx].sum() + S[node.SIdx],L[LIdx].set().count()));
    }
  } while (!nextHeap.empty());


//  cout << endl << endl << "Solution" << endl;
//  for (size_t i=0;i<LParam.size();i++) {
//    cout << std::setw(2) << i << ": Sum = " << std::setw(2) << LParam[i].sum() << " Bits: " << getBits(LParam[i].set(),N,N) << endl;
//  }
//
//  cout << endl << endl << "Reverse Solution" << endl;
//  uint64_t sum=0;
//  for (int i=0;i<N;i++) {
//    sum += S[i];
//  }
//  for (size_t i=0;i<LParam.size();i++) {
//    cout << std::setw(2) << i << ": Sum = " << std::setw(2) << LParam[i].reverseSum(sum) << " Bits: " << getBits(LParam[i].set(),N,N) << endl;
//  }

  return LParam.size();
}


} // end namespace

