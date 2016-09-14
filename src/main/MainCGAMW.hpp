/*
 * CIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef CGA_HPP_
#define CGA_HPP_


#include "Main.hpp"
#include "../globals.hpp"
#include "../partition/CGA_MW.hpp"
#include <sstream>

class MainCGA_MW : public Main {

protected :
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {


			return partition::executeCGA_MW(problem, stats);
	}
	SolutionMethodKPartition getSolutionMethod() {
		return CGA_MW;
	}

};

#endif /* CGA_HPP_ */
