/*
 * Horowitz_Sahni_Cardinality.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: ethan
 */

#include "Horowitz_Sahni_Cardinality.hpp"

#include <algorithm>

namespace ss {

	size_t generateSetsHSCardinality(const uint64_t S[MAXN], const int n,
																	 const uint64_t lower, const uint64_t upper,
																	 const size_t MAX_CARDINALITY, vector<SetNodeBitset> &allSets)

	{
		vector<SetNodeBitset> a;          // all subset sums from first half of numbers
		vector<SetNodeBitset> b;          // all subset sums from second half of numbers

		generateAllSets(S, n, a, 0    , n/2-1  , MAX_CARDINALITY, true);   // combinations of 1st half of numbers (sorted ascending)
		generateAllSets(S, n, b, n/2  , n-1    , MAX_CARDINALITY, false);  // combinations of 2nd quarter of numbers (sorted descending)

		size_t bAnchor = 0;                              			// bPtr starts searching from here

		for (size_t aPtr=0; aPtr < a.size(); aPtr++) {
			const uint64_t &aValue = a[aPtr].sum();        			// For readability, the sum of the a subset.
			const size_t aCardinality = a[aPtr].set().count();	// The cardinality of the a subset.

			while (bAnchor < b.size() &&                   			// Move bAnchor to highest value that might fit
						 aValue + b[bAnchor].sum() > upper) {    			// Keep looking at smaller number in b until we
				bAnchor++;                                   			// are under upper bound. (remember b is descending)
			}

			// Now add all aValue + b[bPtr] that are within range
			for (size_t bPtr=bAnchor;(bPtr < b.size()) && (aValue + b[bPtr].sum() >= lower); bPtr++) {

				// Only add if subset cardinality is no greater than MAX_CARDINALITY
				if (aCardinality + b[bPtr].set().count() <= MAX_CARDINALITY) {
					allSets.push_back(SetNodeBitset(aValue + b[bPtr].sum(),
																			a[aPtr].set() | b[bPtr].set()));
				}
			}
		}

		return allSets.size();                                 // return the number of subset sums stored in ALLSUMS
	}


	// ---------------------------------------------------------------------------------------
	// Helper Functions for generating half sets ONLY of cardinality from 0 to MAX_CARDINALITY
	// TODO: For MAX_CARDINALITY, only include the set if it is in the range of lower and upper
	// 			 bound since it can only be included in the final HS solution if it is combined
	//			 with the empty set
	// ---------------------------------------------------------------------------------------

	/*
		 GENSETS takes an array of integers, an array of subsets to fill, the first
		 index and last indices into the integer array, the sum of the numbers
		 included so far, and the characteristic function of the elements included so
		 far. It generates all subsets of the remaining numbers in the NUMbers array
		 from index FIRST to index LAST, storing them in the SET array.  It places the
		 numbers in consecutive locations indexed by the global variable NEXTSUM, leaving
		 in NEXTSUM the number of sets generated. */

	void generateAllSets (const uint64_t S[], vector<SetNodeBitset> &sets,
												int first, int last, uint64_t cursum, DynamicBitset &curSet,
												size_t cardinality, const size_t MAX_CARDINALITY)

	{
		if (first > last) {                                    // set is completed

			sets.push_back(SetNodeBitset(cursum,curSet));

		} else {                                               // set not yet completed

			// exclude next element
			generateAllSets (S, sets, first+1, last, cursum, curSet, cardinality, MAX_CARDINALITY);


			if (cardinality < MAX_CARDINALITY) {	// If less than MAX_CARDINALITY
				// Include next element
				curSet[first] = true;
				generateAllSets (S, sets, first+1, last, cursum + S[first], curSet, cardinality+1, MAX_CARDINALITY);
				curSet[first] = false;
			}
		}
	}

	// Helper function, provides missing variables
	// also sorts sets before returning
	void generateAllSets(const uint64_t S[], const int N, vector<SetNodeBitset> &sets,
											 int first, int last, const size_t MAX_CARDINALITY, bool sortAscending) {
		DynamicBitset curSet(N);                                   // The number of total input elements
		generateAllSets(S,sets,first,last,0ll,curSet,0,MAX_CARDINALITY);
		if (sortAscending) {
			std::sort( sets.begin(), sets.end() ,std::less<SetNodeBitset>());  // Sort in ascending order
		} else {
			std::sort( sets.begin(), sets.end() ,std::greater<SetNodeBitset>());  // Sort in descending order
		}

	}

}
