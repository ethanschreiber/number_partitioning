/*
 * InclusionExclusion.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: ethan
 */

#ifndef INCLUSIONEXCLUSION_HPP_
#define INCLUSIONEXCLUSION_HPP_
#include "../SubsetSum.hpp"

// -------------------
// Inclusion/Exclusion
// -------------------

namespace ss {
	size_t incExc(const uint64_t S[], const int n, const uint64_t sumRemaining,
												const uint64_t lower, const uint64_t upper,
												vector<SetNodeBitset> &allSets);

	size_t generateSetsIESimpleDominance(const uint64_t S[], const int n, const uint64_t sumRemaining,
												const uint64_t lower, const uint64_t upper,
												vector<SetNodeBitset> &allSets);
	size_t generateSetsIEMoffittDominance(const uint64_t S[], const int n, const uint64_t sumRemaining,
												const uint64_t lower, const uint64_t upper,
												vector<SetNodeBitset> &allSets);
}


#endif /* INCLUSIONEXCLUSION_HPP_ */
