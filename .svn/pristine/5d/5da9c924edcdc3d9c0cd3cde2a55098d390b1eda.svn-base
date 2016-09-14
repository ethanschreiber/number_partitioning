/*
 * Horowitz_Sahni.hpp
 *
 *  Created on: Feb 4, 2014
 *      Author: ethan
 */

#ifndef HOROWITZ_SAHNI_CARDINALITY_HPP_
#define HOROWITZ_SAHNI_CARDINALITY_HPP_

#include "../SubsetSum.hpp"
// ------------------
// Horowitz and Sahni
// ------------------

namespace ss {
	size_t generateSetsHSCardinality(const uint64_t S[], const int n,
																	 const uint64_t lower, const uint64_t upper,
																	 const size_t maxCardinality, vector<SetNodeBitset> &allSets);


	// ---------------------------------------------------------------------------------------
	// Helper Functions for generating half sets ONLY of cardinality from 0 to MAX_CARDINALITY
	// ---------------------------------------------------------------------------------------
	void generateAllSets(const uint64_t S[], const int N, vector<SetNodeBitset> &sets,
	             	 	 	 	 int first, int last, const size_t MAX_CARDINALITY, bool sortAscending);

}
#endif /* HOROWITZ_SAHNI_HPP_ */
