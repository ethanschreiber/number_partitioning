/*
 * ciw.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef MAIN_MOFFITT_HPP_
#define MAIN_MOFFITT_HPP_


#include "Main.hpp"
#include "../globals.hpp"
#include "../partition/Moffitt.hpp"
#include <sstream>

class MainMoffitt : public Main {

protected :
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {
		return partition::executeMoffitt(problem,stats);
	}

	SolutionMethodKPartition getSolutionMethod() {
		return MOFFITT;
	}
};

#endif /* MAIN_MOFFITT_HPP_ */
