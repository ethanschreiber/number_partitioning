/*
 * DynamicProgramming.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: ethan
 */

#ifndef DYNAMICPROGRAMMING_HPP_
#define DYNAMICPROGRAMMING_HPP_

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include "../partition/PartitionUtils.hpp"

namespace ss {
typedef boost::dynamic_bitset<uint64_t> DPBitset;
uint64_t executeDPPartition(const partition::PartitionProblem &problem, ProblemStats &stats);


// --------------------------------------
// Converts idx to value and value to idx
// --------------------------------------
//
// Largest value is at idx 0, then 1, 2 all the way to ROW_SIZE-1 which contains the smallest.
// For example, if ROW_SIZE=10 (0..9), and the values are from 3 to 12, then converting value to idx:
//               Value                    IDX
// idxValue(10, 3, 3)  = 10 - 1 - 3  + 3 = 9 (smallest is last)
// idxValue(10, 3, 7)  = 10 - 1 - 7  + 3 = 5
// idxValue(10, 3, 12) = 10 - 1 - 12 + 3 = 0 (largest is first)
//
// And back from idx to value:
//                IDX                  VALUE
// idxValue(10, 3, 9) = 10 - 1 - 9 + 3 = 3 (smallest is last)
// idxValue(10, 3, 5) = 10 - 1 - 5 + 3 = 7
// idxValue(10, 3, 0) = 10 - 1 - 0 + 3 = 12 (largest is first)
inline size_t idxValue(const size_t ROW_SIZE, const uint64_t MIN_VALUE, const uint64_t IDX_VALUE) {
	return ROW_SIZE - 1 - IDX_VALUE + MIN_VALUE;
}

}
#endif /* DYNAMICPROGRAMMING_HPP_ */
