/*
 * Extended_Schroeppel_Shamir.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: ethan
 */

#ifndef EXTENDED_SCHROEPPEL_SHAMIR_HPP_
#define EXTENDED_SCHROEPPEL_SHAMIR_HPP_

#include "../SchroeppelShamir.hpp"

namespace ss {

// ===================================
// Extended Schroeppel and Shamir Sets
// ===================================

// allSets has a vector of ints
size_t ESSSets(const uint64_t S[MAXN], const int n, const uint64_t lower, const uint64_t upper,
                      vector<SetNodeVector> &allSets);


// Given a set of n input numbers S, a lower and upper bound,
// this generates all sets within lower and upper and puts them in allSets.
// This is a template function, ContainerT must have member functions
// push_back(SSSetNode) and and size().
// Furthermore, the type of the data being stored by the template must be
// SSSetNode. i.e. std::vector<SSSetNode> would be a valid type to pass for allSets

template <typename ContainerT>
size_t ESSSets(const uint64_t S[MAXN], const int n,
                      const uint64_t lower, const uint64_t upper,
                      ContainerT &allSets)
{
  vector<SetNodeBitset> a;    // all subset sums from first quarter of numbers
  vector<SetNodeBitset> b;    // all subset sums from second quarter of numbers
  vector<SetNodeBitset> c;    // all subset sums from third quarter of numbers
  vector<SetNodeBitset> d;    // all subset sums from fourth quarter of numbers

  SSMinHeap abheap;    				// Min Heap Comparator
  SSMaxHeap cdheap;   				// Max Heap Comparator

  deque<SSHeapNode> cdlist;   // list of subset sums removed from CD heap

  // Initialize the 4 vectors and 2 heaps
  initialize(S,n,a,b,c,d,abheap,cdheap,upper);

  // --------------------------------------------------------------------------
  // Now fill in the allSets vector
  // --------------------------------------------------------------------------

  while (!abheap.empty()) {                      // until AB heap is empty

    // remove sets with sums above upper bound from cdlist
    while (!cdlist.empty() && abheap.front().sum() + cdlist.front().sum() > upper) {
      cdlist.pop_front();
    }

    // add new sets with sums within bounds to cdlist
    uint64_t topsum = cdheap.front().sum() + abheap.front().sum();     //m_sum of elements on top of heaps
    while (!cdheap.empty() && topsum >= lower) {         // new set within bounds

      if (topsum <= upper) {                             // new set within bounds
        cdlist.push_back(cdheap.front());                // add heap element to list
      }

      if (cdheap.front().y() > 0) {                        // there is another element in this column
        SSHeapNode newElement(cdheap.front().x(),          // keep same row
                              cdheap.front().y() - 1,      // replace with next combination in column
                              c,d);                      // m_sum is c[x] + d[y-1]

        cdheap.pop();                                    // Remove old head
        cdheap.push(newElement);                         // add new element to heap vector

      } else {                                           // no more elements in this column
        cdheap.pop();                                    // Remove old head
      }

      topsum = cdheap.front().sum() + abheap.front().sum();      // m_sum of two top elements of heaps
    }

    // Combine each element from cdlist with the top element from
    // abheap, these are all within range
    for (size_t index=0;index<cdlist.size();index++) {
      uint64_t subsum = cdlist[index].sum() + abheap.front().sum();  // first subset m_sum
      allSets.push_back(SetNodeBitset(subsum, c[cdlist[index].x()].set() |
                                       d[cdlist[index].y()].set() |
                                       a[abheap.front().x()].set() |
                                       b[abheap.front().y()].set()));
    }

    // If there is another element from the b set, add it.
    // Then, remove top element from abheap.
    if (abheap.front().y() < b.size() - 1) {          // there is another element in this column
      SSHeapNode newNode(abheap.front().x(),            // new element being added to heap
                       abheap.front().y() + 1,        // replace with next combination in this column
                       a,b);

      abheap.pop();                            // Pop old value
      abheap.push(newNode);                    // add new element to heap vector
    } else {
      abheap.pop();
    }
  }

  return allSets.size();                                 // return the number of subset sums stored in ALLSUMS
}


// ===================================
// Extended Schroeppel and Shamir Sums
// ===================================

size_t ESSSums(const uint64_t S[], int n,uint64_t lower, uint64_t upper,
                     vector<uint64_t> &allsums);

} // end namespace

#endif /* SCHROEPPEL_SHAMIR_HPP_ */
