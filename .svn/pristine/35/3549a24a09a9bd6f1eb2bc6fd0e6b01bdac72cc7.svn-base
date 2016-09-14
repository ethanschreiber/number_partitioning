/*
 * Schroeppel_Shamir.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "Extended_Schroeppel_Shamir.hpp"
#include  <iostream>
#include <iomanip>
using std::cout;
using std::endl;

namespace ss {

// ===================================
// Extended Schroeppel and Shamir Sets
// ===================================

size_t ESSSets(const uint64_t S[MAXN], const int n,
					 const uint64_t lower, const uint64_t upper,
					 vector<SetNodeVector> &allSets) {
  vector<SetNodeBitset> allSetsBitSet;

  ESSSets(S,n,lower,upper,allSetsBitSet);  // Call Template function
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




// ===================================
// Extended Schroeppel and Shamir Sums
// ===================================


// This uses Schroeppel and Shamir!
//     ALLSETS takes an array NUM of integers sorted in increasing order,
//     the length N of the vector, the SUM of all the numbers, a LOWER
//     bound and an UPPER bound on the subset sums, and an array ALLSUMS.
//     It computes all possible subset sums of the given numbers, and
//     loads them into the array ALLSUMS, in the order they are generated.
//     It returns the number of subset sums.

size_t ESSSums(const uint64_t S[MAXN], const int n, const uint64_t lower, const uint64_t upper,
             	 	 	 	 	 vector<uint64_t> &allsums)

{
  vector<uint64_t> a;       // all subset sums from first quarter of numbers
  vector<uint64_t> b;       // all subset sums from second quarter of numbers
  vector<uint64_t> c;       // all subset sums from third quarter of numbers
  vector<uint64_t> d;       // all subset sums from fourth quarter of numbers

  SSMinHeap abheap;    			// Min Heap Comparator
  SSMaxHeap cdheap;   			// Max Heap Comparator

  deque<uint64_t> cdlist; 	// list of subset sums removed from CD heap


  // Initialize the 4 vectors and 2 heaps
  initialize(S,n,a,b,c,d,abheap,cdheap,upper);

  // --------------------------------------------------------------------------
  // Now fill in the allsums vector
  // --------------------------------------------------------------------------

  while (!abheap.empty()) {                      // until AB heap is empty

    // first remove sets with sums above upper bound from list
    while (!cdlist.empty() && abheap.front().sum() + cdlist.front() > upper) {
      cdlist.pop_front();
    }
    // add new sets to list with sums within bounds
    uint64_t topsum = cdheap.front().sum() + abheap.front().sum();     //sum of elements on top of heaps
    while (!cdheap.empty() && topsum >= lower) {  // new set within bounds
      if (topsum <= upper) {                      // new set within bounds
        cdlist.push_back(cdheap.front().sum());     // add heap element to list
      }

      if (cdheap.front().y() > 0) {                        // there is another element in this column
        SSHeapNode newElement(cdheap.front().x(),            // keep same row
                            cdheap.front().y() - 1,        // replace with next combination in column
                            c,d);                        // sum is c[x] + d[y]

        cdheap.pop();                                              // Remove old head

        cdheap.push(newElement);                                   // add new element to heap vector

      } else {                                                          // no more elements in this column
        cdheap.pop();                                              // Remove old head
      }

      topsum = cdheap.front().sum() + abheap.front().sum();      // sum of two top elements of heaps
    }


    for (size_t index=0;index<cdlist.size();index++) {
      uint64_t subsum = cdlist[index] + abheap.front().sum();  // first subset sum
      if (subsum >= lower && subsum <= upper) {     // subset sum within bounds
        allsums.push_back(subsum);                  // store subset sum in result vector
      }
    }

    if (abheap.front().y() < b.size() - 1) {            // there is another element in this column
      SSHeapNode newHeapNode(abheap.front().x(),               // replace with next combination in this column
                             abheap.front().y() + 1,
                             a,b);

      abheap.pop();
      abheap.push(newHeapNode);
    } else {
      abheap.pop();
    }
  }
  return allsums.size();                                 // return the number of subset sums stored in ALLSUMS
}

} // end namespace
