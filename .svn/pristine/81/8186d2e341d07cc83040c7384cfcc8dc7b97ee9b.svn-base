/*
 * BinPackingUtils.hpp
 *
 *  Created on: Jul 17, 2012
 *      Author: ethan
 */

#ifndef BINPACKINGUTILS_HPP_
#define BINPACKINGUTILS_HPP_

#define MAXN 70             /* maximum number of elements to be partitioned */
#include "../globals.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>  // for std::accumulate
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <stdexcept>


using std::cout;
using std::endl;
using std::string;
using std::vector;


// ---------------------
// PackingOptions Struct
// ---------------------
struct PackingOptions {

	string inputFilename;
	int solutionMethodInt;
  bool printSolution;
  bool useSchroeppelShamir;
  bool useSchroeppelShamirBinarySearch;
  bool useLDS;
  bool isVerbose;
  bool removePairs;
  bool classicSort; // Should we use the classic korf sort?
  int bufferSize;   // The bin completion buffer size (only if !isClassicSearch && !useSchroeppelShamir

  PackingOptions() :
  	inputFilename(UNSET_STRING),
  	solutionMethodInt(UNSET_INT),
  	printSolution(false),
  	useSchroeppelShamir(false),
  	useSchroeppelShamirBinarySearch(false),
  	useLDS(false),
  	isVerbose(false),
  	removePairs(false),
  	classicSort(false),
  	bufferSize(50){

  }
};


// ===============
// Packing Problem
// ===============
struct PackingProblem {
  int N;                    // The number of elements.
  uint64_t *S;               // The array of N elements.
  uint64_t sum;              // The sum of the elements.
  string problemName;

  PackingProblem(int N0, uint64_t S0[], uint64_t sum0, string problemName0 );
  PackingProblem(const string &filename);
  PackingProblem(const PackingProblem &problem);
  PackingProblem() {}
  void readElements(std::ifstream &inFile);
  virtual ~PackingProblem();
  void load(string filename);
};

// ===================================================================
// A structure to contain the information defining an instance
// of a bin packing problem.  The struct also contains functions to
// read the data from a specification file.
// ===================================================================
struct BinPackingProblem : public PackingProblem {
  uint64_t capacity;         // The capacity of a bin.

  int lowerBound;						// Optional: The lower bound for the problem

  // --------------------------------------------------------------------------
  // Constructor reads the file from standard input. A specification for
  // the file is located here:
  //    http://scip.zib.de/doc/examples/Binpacking/READER.html
  // The specification is as follows:
  //  line 1   : [Problem name]
  //  line 2   : [Capacity] [# of items] [Value of known feasible solution]
  //  line 3 on: [Size of each element one per line]
  // --------------------------------------------------------------------------

  BinPackingProblem(const uint64_t capacity, int lowerBound=-1);  // use this along with reset to create BinPackingProblem
  BinPackingProblem(const string &filename, int lowerBound=-1);
  BinPackingProblem(const PackingProblem &partitioningProblem, uint64_t capacity, int lowerBound=-1);

  void load(string filename);
  void reset(uint64_t newS[MAXN], size_t newN);


  string toString() const;
};



// *****************************
// Bin Packing Utility Functions
// *****************************

// ============================================================================
// BFD takes the array of numbers, and the bins, and puts the next item in the
// first bin it fits into.  The bins are kept sorted in decreasing order of
// size.  It returns the total number of bins used, and initializes SOLUTION to
// the bfd solution.
// ============================================================================

uint64_t BFD(uint64_t items[MAXN],int solution[MAXN],int numItems,uint64_t capacity);

// no solution array
uint64_t BFD(uint64_t items[MAXN],int numItems,uint64_t capacity);

// ============================================================================
// L2 takes an Array of elements, and the index of the FIRST and LAST elements,
//  and returns a lower bound on the number of bins needed for this problem,
//  based on the wasted-space calculation.
// ============================================================================
uint64_t L2(const uint64_t S[MAXN],         // The elements to be assigned
       int N,               // The number of elements
       uint64_t capacity,        // The capacity of a bin
       uint64_t sum);            // Sum of all of the elements.



// ============================================================================
// SetNodeVector Class - This is the SetNode class that stores elements as a
// vector of uint64_t. Not to be confused with SetNodeBitset.
// ============================================================================
class SetNodeVector {              // type of a set of elements
public :
  static const uint64_t SENTINEL_VALUE = -1;   // For ending loops, might be last element

 uint64_t sum;                     // the sum of the elements in the set
 vector<uint64_t> member;          // The elements, i.e., their sizes

 SetNodeVector(uint64_t sum0, int numIncluded, uint64_t member0[]) : sum(sum0), member(numIncluded) {

   for (int i = 0; i < numIncluded; i++) {        // for each element of set
     member[i] = member0[i];                      // store element in set
   }
 }
 SetNodeVector(int numIncluded, uint64_t member0[]) : sum(0), member(numIncluded) {

   for (int i = 0; i < numIncluded; i++) {        // for each element of set
     member[i] = member0[i];                      // store element in set
     sum+=member0[i];
   }
 }
 SetNodeVector() : sum(0),member() {
 }


 const vector<uint64_t> &elements() const {
   return member;
 }

 uint64_t getSum() const {
   return sum;
 }

 // Puts firstItem in the list and resizes the list
 // to size 1
 void replaceWith(uint64_t firstItem) {
   member.resize(1,firstItem);
   member[0] = firstItem;
 }

 void push_back(uint64_t item) {

   sum += item;
   member.push_back(item);
 }

 void push_sentinel() {
   member.push_back((uint64_t) SENTINEL_VALUE);
 }

 void pop_back() {
   if (member.back() != SENTINEL_VALUE) {
     sum -= member.back();
   }
   member.pop_back();
 }

 void pop_back(int numToPop) {
   for (int i=0;i<numToPop && !empty();i++) {
     pop_back();
   }
 }

 uint64_t back() const {
   return member.back();
 }

 bool empty() const {
   return member.empty();
 }

 void resize(size_t size) {
   member.resize(size);
 }

 int cardinality() const {
   if ((!member.empty()) && member[member.size()-1] == SENTINEL_VALUE) {
     return member.size() - 1;
   } else {
     return member.size();
   }
 }

 const uint64_t& operator[](size_t i) const{
   return member[i];
 }

 void clear() {
   member.clear();
   sum = 0;
 }

 string toString() const {
   std::ostringstream out;
   out << "[card = " << cardinality() << " members = {";
   for (int i=0;i<cardinality();i++) {
     if (i != 0) {
       out << ",";
     }
     out << member[i];
   }
   out << "} sum = " << sum << "]";
   return out.str();
 }
};

struct SetNodeCardinalityComparatorClassic {

  //   This is the original one
    bool operator() (const SetNodeVector &node1, const SetNodeVector &node2) {
      return node1.sum > node2.sum;
    }

};

struct SetNodeCardinalityComparator {

  // is node1 < node2
  bool operator() (const SetNodeVector &node1, const SetNodeVector &node2) {
    if (node1.sum == node2.sum) {                         // If same sum
      int card1 = node1.cardinality();
      int card2 = node2.cardinality();
      if (card1 == card2) {                               // If same cardinality
        for (int i=card1-1;i>=0;i--) {
          if (node1.member[i] != node2.member[i] ) {
            return node1.member[i] < node2.member[i];
          }
        }
        return false;     // They are equal, so not <
      } else {                                            // If not same cardinality, smaller cardinality comes first
        return card1 < card2;
      }
    } else {                                              // If not same sum, largest sum comes first
      return node1.sum > node2.sum;
    }
  }
};

// ============================================================================
// SumVector Class - This is a static array which maintains the sum of its elements
//                   as well. T must implement the "+" operator.
// ============================================================================
template <class T> class SumArray  {
private :
  T      m_array[MAXN]; // Underlying array for storing data
  T     *m_begin;       // A pointer to the first element we care about
  T      m_sum;         // The sum of the elements in the array
  size_t m_size;        // The number of elements (starting at m_begin)

public :
  SumArray() : m_begin(m_array), m_sum(0), m_size(0) {
  }

  inline T* begin() const {
    return m_begin;
  }

  inline T lastElement() const {
    return m_array[m_size-1];
  }

  inline size_t size() const {
    return m_size;
  }

  inline bool empty() const {
    return m_size == 0;
  }

  inline T sum() const {
    return m_sum;
  }

  void push_back(T val) {
    m_array[m_size] = val;
    m_size++;
    m_sum += val;
  }

  T pop_front() {
    T val = *m_begin;
    m_sum -= val;  // Subtract first element from sum
    m_begin += 1;       // Have begin point to next element
    m_size--;           // 1 fewer element now

    return val;
  }

  T& operator[](size_t index) {
      return *(m_begin+index);
  }
};

#endif /* BINPACKINGUTILS_HPP_ */




