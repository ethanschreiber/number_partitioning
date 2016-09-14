/*
 * CGA.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: ethan
 */

#ifndef CGA_HPP_
#define CGA_HPP_
#include <stdint.h>
// All functions assume S or problem.S are sorted in descending order

namespace ss {
  // Exact
  uint64_t executeCGA2(const uint64_t S[], const int N, const uint64_t sum);

  // Heuristic
  uint64_t executeGreedy(const uint64_t S[], const int N, const int K);
}
#endif /* CGA_HPP_ */
