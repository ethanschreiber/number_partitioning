/*
 * CIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef BSBCP_HPP_
#define BSBCP_HPP_

#include "MainBS.hpp"
#include <sstream>
#include "../partition/BinarySearchBCP.hpp"


class MainBSBC : public MainBS {

protected :

	uint64_t executeBSPartition(BinPackingProblem &binPackingProblem,	uint64_t minCapacity, uint64_t maxCapacity,
			int K,uint64_t maxElement,	const PackingOptions &packingOptions) {

		binary_search::BinarySearchBCP bsbcp;
		return bsbcp.search(binPackingProblem,minCapacity,maxCapacity+1,
											 K,maxElement,packingOptions);
	}
	SolutionMethodKPartition getSolutionMethod() {
		return BSBCP;
	}



};

#endif /* BSBCP_HPP_ */

