/*
 * Rich.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: ethan
 */

#ifndef RICH_HPP_
#define RICH_HPP_
#include "PartitionUtils.hpp"

#define MAXK 13              /* maximum number of sets being partitioned into */

namespace partition {
	uint64_t executeSNP(const partition::PartitionProblem &problem, ProblemStats &stats);
}




#endif /* RICH_HPP_ */
