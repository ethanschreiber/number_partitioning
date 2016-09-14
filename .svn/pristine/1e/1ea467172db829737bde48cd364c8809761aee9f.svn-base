/*
 * PartitionUtils.hpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ethan
 */

#ifndef PARTITIONUTILS_HPP_
#define PARTITIONUTILS_HPP_

#include "../pack/PackingUtils.hpp"
#include "../ss/SubsetSum.hpp"
#include "../ss/KK.hpp"
#include "../ss/CGA.hpp"

namespace partition {

typedef vector<ss::SetNodeBitset> SetVector;
typedef deque<ss::SetNodeBitset> SetDeque;
typedef SetVector::const_iterator SetVectorIt;
typedef SetDeque::const_iterator SetDequeIt;

typedef std::pair<SetVectorIt,SetVectorIt> SetIteratorPair;

// ---------------------
// PartitionOptions Struct
// ---------------------

struct PartitionOptions {

  string          inputFilename;
  string					solutionFilename;		// This is for extended subset sum, find out solution value to get upper bound
  int             solutionMethodInt;
  int             inputK;
  int             problemIdx;
  uint64_t        minCapacity;
  uint64_t        maxCapacity;
  bool            isLowCardinality;		// For CIW, if true, uses low cardinality version
  size_t          numSets;


  PartitionOptions() :
  	inputFilename     (UNSET_STRING),
  	solutionFilename  (UNSET_STRING),
    solutionMethodInt (UNSET_INT),
    inputK            (UNSET_INT),
    problemIdx        (UNSET_INT),
    minCapacity       (UNSET_UINT64_T),
    maxCapacity       (UNSET_UINT64_T),
  	isLowCardinality  (false),
  	numSets           ((isLowCardinality) ? 1000 : 20000)
  {

  }
};


// ==========================
// PartitioningProblem Struct
// ==========================
struct PartitionProblem : public PackingProblem {

  int K;    // The number of partitions

  // --------------------------------------------------------------------------
  // A Partitioning problem specification is similar to the BinPackingProblem
  // described below. The specification is as follows:
  //
  // line 1: [Problem name]
  // line 2a: [# of items] [sum]    (k not specified, specified outside of file)
  // line 2b: [k] [# of items] [k]  (k is specified)
  // line 3 on: [Size of each element one per line]
  // --------------------------------------------------------------------------

  PartitionProblem(int N0, int K0, uint64_t S0[], uint64_t sum0, string problemName0 );
  PartitionProblem(const string &filename);
  PartitionProblem(const PartitionProblem &problem);
  PartitionProblem(const PackingProblem &problem, int K0);
  void load(const string &filename);
  string toString() const;

};

	//      To round up q = x/y :
	//      - q = (x + y - 1) / y;
	//      - or (avoiding overflow in x+y)
	//        q = 1 + ((x - 1) / y); // if x != 0
template <class T>
T divCeiling (T x, T y) {
	 return 1 + ((x - 1) / y);
	}


// ================
// Utility Function
// ================
template <typename BitsetType>
string getBits(const BitsetType &bitset,int groupSize,int numBits) {
  std::ostringstream out;
  for (int i=numBits-1;i>=0;i--) {
    out << bitset[i];
    if (i % groupSize == 0 && i > 0) {
      out << " ";
    }
  }
  return out.str();
}

template <typename BitsetType>
string getBitsIndices(const BitsetType &bitset,int numBits) {
  std::ostringstream out;

  for (int i=0;i<numBits;i++) {
    if (bitset[i]) {
    	out << std::setw(2) << i << " ";
    }
  }
  return out.str();
}

	inline uint64_t getLowerBound(const uint64_t sum, const int K) {
		return divCeiling(sum, (uint64_t) K);  // Best we can do is divide evenly
	}

	inline uint64_t getLowerBound(const PartitionProblem &problem) {
		return getLowerBound(problem.sum,problem.K);
	}

	// Used to compute a new CMin whenever ub changes.
	inline uint64_t computeCMin(uint64_t elementsSum, int K, uint64_t ub) {
	  uint64_t val =  (K-1) * (ub-1);
	  if (elementsSum > val) {
	    return elementsSum - val;
	  } else {
	    return 0;
	  }
	}

	inline uint64_t getUpperBound(const PartitionProblem &problem) {
	 uint64_t kkVal = kk(problem.S,problem.N,problem.K,problem.sum);
	 return kkVal;
	}

	inline string spaces(int K) {
		std::ostringstream out;
		for (int i=0;i<K;i++) {
			out << " ";
		}
		return out.str();
	}


} 	// end namespace


#endif /* PARTITIONUTILS_HPP_ */
