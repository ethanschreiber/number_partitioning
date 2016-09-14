/*
 * CIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef BSBC_HPP_
#define BSBC_HPP_

#include "MainBS.hpp"
#include <sstream>
#include "../partition/BinarySearchBC.hpp"

class MainBSBC : public MainBS {

protected :
	uint64_t executeBSPartition(BinPackingProblem &binPackingProblem,	uint64_t minCapacity, uint64_t maxCapacity,
			int K,uint64_t maxElement,	const PackingOptions &packingOptions) {

		binary_search::BinarySearchBC bsbc;
		return bsbc.search(binPackingProblem,minCapacity,maxCapacity+1,
											 K,maxElement,packingOptions);
	}

	SolutionMethodKPartition getSolutionMethod() {
		return BSBC;
	}
};

#endif /* BSBC_HPP_ */

