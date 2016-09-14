/*
 * CIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef CIW_HPP_
#define CIW_HPP_


#include "Main.hpp"
#include "../globals.hpp"
#include "../partition/PartitionIterativeWeakening.hpp"
#include "../partition/PartitionIterativeWeakeningLowCardinality.hpp"
#include <sstream>

class MainCIW : public Main {
private:
	bool m_isLowCardinality;
protected :
	uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
			const PackingOptions &packingOptions, ProblemStats &stats) {


		if (m_isLowCardinality) {
			return partition::executeCIWLowCardinality(problem,stats,partitionOptions.numSets);
		} else {
			return partition::executeCIW(problem, stats,partitionOptions.numSets);
		}
	}

	void init(partition::PartitionOptions &partitionOptions) {
		m_isLowCardinality = partitionOptions.isLowCardinality;
	}
	SolutionMethodKPartition getSolutionMethod() {
		if (m_isLowCardinality) {
			return CIW_LC;
		} else {
			return CIW;
		}
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

#endif /* CIW_HPP_ */
