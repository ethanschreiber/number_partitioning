/*
 * globals.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: ethan
 */

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_
#include <limits.h>
#include <float.h>
#include <string>
#include <stdint.h>
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

// In general to mark an unset value
const int       UNSET_INT       = INT_MAX;
const string    UNSET_STRING    = "";
const uint32_t  UNSET_UINT32_T  = UINT_MAX;
const uint64_t  UNSET_UINT64_T  = ULONG_MAX;
const double    UNSET_DOUBLE    = DBL_MAX;

// const int       UNSET_INT       = std::numeric_limits<int>::max();
// const string    UNSET_STRING    = "";
// const uint32_t  UNSET_UINT32_T  = std::numeric_limits<uint32_t>::max();
// const uint64_t  UNSET_UINT64_T  = std::numeric_limits<uint64_t>::max();
// const double    UNSET_DOUBLE    = std::numeric_limits<double>::max();


static string methodsToString(const string SUFFIXES[], const size_t NUM_METHODS) {
	std::ostringstream out;
	for (size_t i=0;i<NUM_METHODS;i++) {
		out << ((i != 0) ? ", " : "")
				<< i << "=" << (SUFFIXES[i].substr(1));	// Remove "_"
	}
	return out.str();
}

// ==============================
// For Multiway Partition Problem
// ==============================
const size_t NUM_K_PARTITION_METHODS = 9;
enum SolutionMethodKPartition {
	RNP 		= 0,	// Recursive Number Partitioning (Korf)
	BSBC		= 1, 	// Binary-Search Bin Completion (Schreiber and Korf)
	BSBCP		= 2,	// Binary-Search Branch-and-Cut-and-Price (Belov)
	MOFFITT	= 3, 	// Moffitt's algorithm
	SNP			= 4, 	// Sequential Number Partitioning (Korf and Schreiber)
	CIW			= 5, 	// Cached Iterative Weakening	(Schreiber and Korf)
	CIW_LC	= 6,	// Cached Iterative Weakening Low Cardinality	(Schreiber and Korf)
	RNP_2009= 7,  // Recursive Number Partitioning (Korf)
	CGA_MW=8		// Complete greedy algorithm multiway
};

const string K_PARTITION_SUFFIXES[NUM_K_PARTITION_METHODS] =
{"_rnp", "_bsbc", "_bsbcp", "_moffitt", "_snp", "_ciw", "_ciwlc", "_rnp2009", "_cga_mw"};

// Prints a list of the Partition Methods
static string kPartitionMethodsToString() {
	return methodsToString(K_PARTITION_SUFFIXES,NUM_K_PARTITION_METHODS);
}

// Get SolutionMethodKPartition from solutionMethodIdx
static inline SolutionMethodKPartition getSolutionMethodKPartition(int solutionMethodIdx) {
	if (solutionMethodIdx < 0 || solutionMethodIdx >= (int) NUM_K_PARTITION_METHODS) {
    cout << "ERROR: Must enter " << kPartitionMethodsToString() << "." << endl;
    exit(0);
	}
	return (SolutionMethodKPartition) solutionMethodIdx;
}

// Get Suffix from method
static inline string getSuffix(SolutionMethodKPartition method) {
	return K_PARTITION_SUFFIXES[(int) method];
}

// ===========================
// For 2-way Partition Problem
// ============================
const size_t NUM_PARTITION_METHODS = 5;
enum SolutionMethodPartition {
	CGA 		= 0,	// Complete Greedy Algorithm (Korf)
	CKK 		= 1, 	// Complete Karmarkar-Karp(Korf)
	HS      = 2,	// Horowitz and Sahni
	SS    	= 3, 	// Schroeppel and Shamir
	DP			= 4, 	// Dynamic Programming
};

const string PARTITION_SUFFIXES[NUM_PARTITION_METHODS] =
{"_cga", "_ckk", "_hs", "_ss", "_dp"};

// Prints a list of the Partition Methods
static string partitionMethodsToString() {
	return methodsToString(PARTITION_SUFFIXES,NUM_PARTITION_METHODS);
}

// Get SolutionMethodPartition from solutionMethodIdx
static inline SolutionMethodPartition getSolutionMethodPartition(int solutionMethodIdx) {
	if (solutionMethodIdx < 0 || solutionMethodIdx >= (int) NUM_PARTITION_METHODS) {
    cout << "ERROR: Must enter " << partitionMethodsToString() << "." << endl;
    exit(0);
	}
	return (SolutionMethodPartition) solutionMethodIdx;
}

// Get Suffix from method
static inline string getSuffix(SolutionMethodPartition method) {
	return PARTITION_SUFFIXES[(int) method];
}

// ==============================
// For Extended Partition Problem
// ==============================
const size_t NUM_EXTENDED_PARTITION_METHODS = 2;
enum SolutionMethodExtendedPartition {
	IE      = 0,	// Inclusion-Exclusion
	ESS    	= 1 	// Extended Schroeppel and Shamir
};

const string EXTENDED_PARTITION_SUFFIXES[NUM_EXTENDED_PARTITION_METHODS] =
{"_ie", "_ess"};

// Prints a list of the Partition Methods
static string extendedPartitionMethodsToString() {
	return methodsToString(EXTENDED_PARTITION_SUFFIXES,NUM_EXTENDED_PARTITION_METHODS);
}

// Get SolutionMethodPartition from solutionMethodIdx
static inline SolutionMethodExtendedPartition getSolutionMethodExtendedPartition(int solutionMethodIdx) {
	if (solutionMethodIdx < 0 || solutionMethodIdx >= (int) NUM_EXTENDED_PARTITION_METHODS) {
    cout << "ERROR: Must enter " << extendedPartitionMethodsToString() << "." << endl;
    exit(0);
	}
	return (SolutionMethodExtendedPartition) solutionMethodIdx;
}

// Get Suffix from method
static inline string getSuffix(SolutionMethodExtendedPartition method) {
	return EXTENDED_PARTITION_SUFFIXES[(int) method];
}




// ==============================
// For Bin Packing Problem
// ==============================
//const size_t NUM_BIN_PACKING_METHODS = 2;
//enum SolutionMethodKPartition {
//	BC		  = 1, 	// Bin Completion (Korf)
//	BCP			= 2,	// Branch and Cut and Price (Belov)
//};
//
//const string K_PARTITION_SUFFIXES[NUM_K_PARTITION_METHODS] =
//{"_rnp", "_bsbc", "_bsbcp", "_moffitt", "_snp", "_ciw", "_ciwlc"};
//
//string kPartitionMethodToString() {
//	std::ostringstream out;
//	for (int i=0;i<NUM_K_PARTITION_METHODS;i++) {
//		out << (i != 0) ? ", " : ""
//				<< i << "=" << K_PARTITION_SUFFIXES[i];
//	}
//	return out.str();
//}
//
//void assertKPartitionIdx(int solutionMethodIdx) {
//	if (solutionMethodIdx < 0 || solutionMethodIdx >= NUM_K_PARTITION_METHODS) {
//    cout << "ERROR: Must enter " << kPartitionMethodToString() << "." << endl;
//	}
//	exit(0);
//}
//
//// Get SolutionMethodKPartition from solutionMethodIdx
//static inline SolutionMethodKPartition getSolutionMethodKPartition(int solutionMethodIdx) {
//	assertKPartitionIdx(solutionMethodIdx);
//	return (SolutionMethodKPartition) solutionMethodIdx;
//}
//
//// Get Suffix from solutionMethodIdx
//static inline string getSuffixKPartition(int solutionMethodIdx) {
//	assertKPartitionIdx(solutionMethodIdx);
//	return K_PARTITION_SUFFIXES[solutionMethodIdx];
//}

// For discerning options
const string SS_STRING        			= "_ss";
const string CLASSIC_SORT_STRING   	= "_oldsort";
const string PAIRS_STRING     			= "_pairs";
const string LDS_STRING       			= "_lds";


// ============================================
// For reporting stats after executing programs
// ============================================
struct ProblemStats {
  string problemName;     // The problem name
  uint64_t numNodes;      // The # of nodes explored to solve the problem
  double time;            // The time used to solve the problem
  uint64_t lowerbound;    // The initial lower bound
  uint64_t sum;           // The sum of the input elements
  int numBins;            // The number of bins used for the solution
  uint64_t maxCapUsed;    // The max capacity used in any bin in the optimal solution
  int numPerfectPairs;    // The number of pairs equaling the bin capacity

  // This is for cached cardinality
  size_t firstCount;			// The number of first bins we look at
  double ssTime;					// The amount of time used by schroeppel and shamir
  int ssCalls;						// The number of ss calls

  double residentMemory;		// The amount of resident memory used for CIE Trees;
  double residentMemorySS;	// The amount of resident memory used for SS;

};

struct CachedPartitioningStats {
  uint64_t numNodes;      // The # of nodes explored to solve the problem
  double timeSS;					// The schroeppel and Shamir time
  size_t firstSetsCount;	// The number of first sets
};

#endif /* GLOBALS_HPP_ */


