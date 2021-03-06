/*
 * BinarySearch.hpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ethan
 */

#ifndef BINARYSEARCH_HPP_
#define BINARYSEARCH_HPP_

#include "../Utils.hpp"
#include "../pack/completions/SSCompletionGenerator.hpp"
#include "../pack/completions/IECompletionGenerator.hpp"

void executeBP( const BinPackingProblem &problem,ProblemStats &stats, const PackingOptions &packingOptions);
void executeBelovBCP( const BinPackingProblem &problem,ProblemStats &stats);

namespace binary_search {

class BinarySearch {

protected :
	virtual void execute( const BinPackingProblem &problem,ProblemStats &stats, const PackingOptions &packingOptions) = 0;
public :
	uint64_t search(BinPackingProblem &problem, uint64_t minCapacity,uint64_t maxCapacity,
	                int numPartitions, uint64_t maxElement,
	                const PackingOptions &packingOptions);
	~BinarySearch() {}

};

} // end namespace


#endif /* BINARYSEARCH_HPP_ */
