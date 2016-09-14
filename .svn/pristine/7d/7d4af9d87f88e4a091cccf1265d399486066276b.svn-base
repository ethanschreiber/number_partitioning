/*
 * CIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef BSBS_HPP_
#define BSBS_HPP_

#include "Main.hpp"
#include <sstream>
#include <math.h>       /* ceil */

#include "../partition/BinarySearchBC.hpp"

// Binary Search Bin Packing (BSBP)includes:
// Binary Search Bin Completetion (BSBC)
// Binary Search Branch-and-Cut-and-Price (BSBCP)

class MainBS : public Main {

protected :

	virtual uint64_t executeBSPartition(BinPackingProblem &binPackingProblem,	uint64_t minCapacity, uint64_t maxCapacity,
			int K, uint64_t maxElement,	const PackingOptions &packingOptions) = 0;

//	binary_search::BinarySearchBC bsbc;
//		return bsbc.search(binPackingProblem,minCapacity,maxCapacity+1,
//							  			 problem.K,maxElement,packingOptions);
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {

		std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order

		// Find the max input element
		uint64_t maxElement = *std::max_element(problem.S,problem.S + problem.N);

		uint64_t minCapacity = partitionOptions.minCapacity;
		uint64_t maxCapacity = partitionOptions.maxCapacity;
		// If the minCapacity was not passed from the command line
		if (minCapacity == UNSET_UINT64_T) {
			// Can't be smaller than ([the sum of the elements] / [# partitions]) and also can't be smaller than maxElement
			minCapacity = std::max((uint64_t) ceil((double) problem.sum / (double) problem.K),maxElement);

			// Also can't be smaller than the kth + k+1st element
			minCapacity = std::max(minCapacity,problem.S[problem.K] + problem.S[problem.K]);
		}

		if (maxCapacity == UNSET_UINT64_T) {
			maxCapacity = kk(problem.S,problem.N,problem.K,problem.sum);
		}

		BinPackingProblem binPackingProblem(problem,maxCapacity);

		return executeBSPartition(binPackingProblem,minCapacity,maxCapacity+1,
											 problem.K,maxElement,packingOptions);
	}



	string processResults(partition::PartitionOptions &partitionOptions, const ProblemStats stats,
			partition::PartitionProblem &problem, const double timeElapsed, const uint64_t result, const int idx) {

		// Set the number of sets to generate to the largest needed so far
		if (idx == 0) {
    	partitionOptions.numSets = stats.firstCount;
    } else{
    	partitionOptions.numSets = std::max(partitionOptions.numSets,stats.firstCount);
    }
		std::ostringstream out;
		out << problem.problemName <<  " "
				<< timeElapsed  <<  " "
				<< result << " "
				<< stats.firstCount << " "
				<< stats.ssCalls << " "
				<< stats.ssTime << " "
				<< std::fixed << std::setprecision(2) << stats.residentMemory << endl;

		return out.str();
	}


};

#endif /* BSBC_HPP_ */

