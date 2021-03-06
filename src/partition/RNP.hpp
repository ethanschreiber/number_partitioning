/*
* RNP.hpp
 *
 *  Created on: Mar 11, 2014
 *      Author: ethan
 */

#ifndef RNP_HPP_
#define RNP_HPP_


#include "../partition/PartitionUtils.hpp"
#include <stdio.h>                                    /* standard I/O library */
#include <math.h>                                     /* mathematical library */
#include <stdint.h>                                   /* for int64_t */
#include <stdlib.h>                                   /* For exit */
#include <iostream>

namespace partition {
	uint64_t executeRNP    (const partition::PartitionProblem &problem, ProblemStats &stats);
	uint64_t executeRNP2009(const partition::PartitionProblem &problem, ProblemStats &stats);

}


#endif /* RNP_HPP_ */

