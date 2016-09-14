/*
 * BinarySearch.hpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ethan
 */

#ifndef BINARYSEARCH_BCP_HPP_
#define BINARYSEARCH_BCP_HPP_

#include "BinarySearch.hpp"

void executeBelovBCP( const BinPackingProblem &problem,ProblemStats &stats);

namespace binary_search {

	class BinarySearchBCP : public BinarySearch {

	protected :
		void execute( const BinPackingProblem &problem,ProblemStats &stats, const PackingOptions &packingOptions);

	};

} // end namespace


#endif /* BINARYSEARCH_HPP_ */
