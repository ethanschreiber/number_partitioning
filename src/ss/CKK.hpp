/*
 * SchroeppelShamirCompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef CKK_HPP_
#define CKK_HPP_

#include "../pack/PackingUtils.hpp"
namespace ss {

// Returns the difference
uint64_t ckk2(uint64_t S[],                       		// array of numbers
            const int n,                            // number of elements in array
            const uint64_t sumRemaining,            // sum of all numbers
            uint64_t best=( ((uint64_t) 1) << 60)); // Don't set this when called externally

uint64_t executeCKK2(uint64_t S[],        // array of numbers
                    const int n,            // number of elements in array
                    const uint64_t sum);    // sum of all numbers
} /* namespace ss */
#endif /* CKK_HPP_ */
