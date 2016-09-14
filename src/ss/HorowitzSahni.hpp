/*
 * HorowitzSahni.hpp
 *
 *  Created on: May 8, 2014
 *      Author: ethan
 */

#ifndef HOROWITZSAHNI_HPP_
#define HOROWITZSAHNI_HPP_
#include "SubsetSum.hpp"
#include "../partition/PartitionUtils.hpp"



using std::vector;
namespace ss {

/**
 * Given a new candidate solution, the sum of S and perfect, return the new best.
 * Returns perfect if the newValue was perfect or 1 less than perfect.
 * If value > perfect, returns value.
 * If value < perfect, returns sum - value
 */
inline uint64_t getPartitionBest(const uint64_t value, const uint64_t perfect, const uint64_t sum) {

  if ((value == perfect) ||
      (value == perfect-1)) {      // if value is perfect, return perfect
    return perfect;
  } else if (value > perfect) {    // If > perfect,
    return value;                  // then this is best
  } else {                         // If < perfect
    return sum - value;            // best is complement of this value
  }
}

/**
 * Given the best solution found so far, and the sum of the input elements S,
 * returns the lower bound
 */
inline uint64_t getLowerBound(const uint64_t best, const uint64_t sum) {
	return sum - (best - 1);
}
/**
 * Run Horowitz and Sahni to solve the (two-way) Partition problem.
 */
uint64_t executeHS2(const partition::PartitionProblem &problem, ProblemStats &stats);

}
#endif /* HOROWITZSAHNI_HPP_ */
