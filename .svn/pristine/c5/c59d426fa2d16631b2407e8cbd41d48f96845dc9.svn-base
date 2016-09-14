/*
 * ciw.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef MAIN_SNP_HPP_
#define MAIN_SNP_HPP_


#include "Main.hpp"
#include "../globals.hpp"
#include "../partition/SNP.hpp"
#include <sstream>

class MainSNP : public Main {

protected :
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {
		return  partition::executeSNP(problem,stats);
	}

	SolutionMethodKPartition getSolutionMethod() {
		return SNP;
	}
};

#endif /* MAIN_SNP_HPP_ */
