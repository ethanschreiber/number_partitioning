/*
 * BinarySearch.hpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ethan
 */

#ifndef BINARYSEARCH_BC_HPP_
#define BINARYSEARCH_BC_HPP_

#include "BinarySearch.hpp"
#include "../pack/BinCompletion.hpp"

namespace binary_search {

class BinarySearchBC : public BinarySearch {

protected :
	void execute( const BinPackingProblem &problem,ProblemStats &stats, const PackingOptions &packingOptions);

};

} // end namespace


#endif /* BINARYSEARCH_HPP_ */
