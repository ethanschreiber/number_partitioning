/*
 * Horowitz_Sahni.hpp
 *
 *  Created on: Feb 4, 2014
 *      Author: ethan
 */

#ifndef EXTENDED_HOROWITZ_SAHNI_HPP_
#define EXTENDED_HOROWITZ_SAHNI_HPP_

#include "../SubsetSum.hpp"
// ------------------
// Horowitz and Sahni
// ------------------
namespace ss {
	size_t EHS(const uint64_t S[], const int n,
												const uint64_t lower, const uint64_t upper,
												vector<SetNodeBitset> &allSets);

	size_t EHS(const uint64_t S[], const int n,
												const uint64_t lower, const uint64_t upper,
												vector<SetNodeVector> &allSets);

}


#endif /* HOROWITZ_SAHNI_HPP_ */
