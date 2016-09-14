/*
 * ciw.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef MAIN_RNP_HPP_
#define MAIN_RNP_HPP_


#include "Main.hpp"
#include "../globals.hpp"
#include "../partition/RNP.hpp"
#include <sstream>

class MainRNP : public Main {

protected :
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {
		return  partition::executeRNP(problem,stats);

	}

	SolutionMethodKPartition getSolutionMethod() {
		return RNP;
	}
};

#endif /* MAIN_RNP_HPP_ */
