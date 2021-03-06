/*
 * CachedIECard.cpp
 *
 *  Created on: Feb 7, 2014
 *      Author: ethan
 */

#include "GenerateCachedIECardinality.hpp"
#include <sstream>
namespace partition {
  // ============================================================================
  //
  // C a c h e d I E C a r d i n a l i t y I E - (CIEC)
  //
  // ============================================================================

  GCIEC::GenerateCachedIECardinality(Partition *partition,
                                     const vector<uint64_t> &S, CachedIETrees *cachedIETrees,
                                     uint64_t elementsSum, int K, uint64_t upperBound, size_t numSets)
  :m_isTmpVerbose(false), m_cachedIETrees(cachedIETrees), m_partition(partition), m_S(S.begin(),S.end()), m_K(K)
  {
  	// -------------------------
  	// Setup debugging variables
  	// -------------------------
    d_nodesExplored.resize(S.size());
    d_timesCalled.resize(S.size());

    for (size_t i=0;i<S.size();i++) {
    	d_nodesExplored[i] = 0;
      d_timesCalled[i] = 0;
    }
  }

  GCIEC::~GenerateCachedIECardinality() {
  }

  // elements has count bits set to one. the elements were shifted by first Idx to the right.
  // count is the number of bits set in elements and targetCount is the number of bits you want set
  //          idx tens:     1100 0000 0000  (this is just 11, 10, 09, ... 01, 00)
  //          idx ones:     1098 7654 3210  (read top to bottom)
  // Example: elements    = 1111 1001 0111
  //          firstIdx    = 8
  //          count       = 9 (number of bits set in elements)
  //          targetCount = 5 (find the idx with 5 bits to the left)
  //
  //          solution = idx + firstIdx = 07 + 8 = 15
  size_t getFirstIdxMax(const ss::DynamicBitset &elements, const size_t firstIdx,
                        const int count,const int targetCount) {

    size_t setBitsRemoved = 0;
    const size_t setBitsToRemove = count - targetCount;	// Have count bits now, want targetCount
    size_t idx = elements.find_first();									// Find idx of first element

//  	cout << "shifted elements: " << getBits(elements,4,elements.size()) << endl;
//  	cout << "first idx       : " << (int) firstIdx << endl
//  			 << "count           : " << count << endl
//  			 << "target count    : " << targetCount << endl
//  			 << "set bits to remo: " << setBitsToRemove << endl << endl;

    while (idx != ss::DynamicBitset::npos &&			// While still bits left, should never happen
           setBitsRemoved < setBitsToRemove) {		// And more to remove
      idx = elements.find_next(idx);							// Find next one, removing first bit
      setBitsRemoved++;														// Count number of bits found
    }

    return firstIdx + idx;		// Offset by firstIdx
  }

  uint64_t GCIEC::generateRestCIE(const CachedIENode *node,					// The current node in the IE tree
                         const int K,									// The idx of the current bin we are filling
                         ss::DynamicBitset &elements,	// Bitset with 1 representing elements still left
                         const uint64_t elementsSum,	// The sum of the elements still left
                         const size_t elementsCount,	// The count of the elements left (# of 1's in elements)
                         const uint64_t pathSum,      // The sum of the elements included so far in this tree recursion
                         size_t cardinality,					// The cardinality of sets in the tree we are searching now
                         size_t firstIdx,							// The first idx we can include (in order to stop duplicate permutations)
                         size_t firstIdxMax,					// The last idx we can include so we can still make enough completions with this cardinality
                         uint64_t lb,									// The minimum sum of a set that could lead to a solution
                         uint64_t partialCost,				// The max sum used so far above this point in the tree
                         uint64_t ub) { 							// The maximum sum of a bin determined by the first set sum (k==m_k)


  	//cout << "  d_nodesExplored(" << getBits(elements,4,elements.size()) + ") ec: " << elementsCount << " K: " << K << "  First idx: " << firstIdx << "  fim: " << firstIdxMax<< endl;

  	d_nodesExplored[cardinality]++;
    // ---------------------------------------------
    // All elements included or excluded (Base Case)
    // ---------------------------------------------
    if (node == NULL) {
//    		cout << "    NODE NULL: Path Sum: " << pathSum << " (lb,ub): " << lb << "   " << ub << endl;

//    	cout << "       genRestCIE lb: " << lb<< endl;
//    	cout << "       genRestCIE ps: " << pathSum << endl;
//    	cout << "       genRestCIE ub: " << ub<< endl << endl;
      if (pathSum >= lb && pathSum < ub) {
//      		cout << "    PARTITION(" << K-1 << ", " << getBits(elements,4,elements.size()) << "...)" << endl;

        ub = m_partition->partition(K-1,elements,elementsSum, elementsCount, cardinality,
                                    firstIdx, std::max(partialCost,pathSum),ub);
      }
    } else {
    	if (cardinality == 1) {
    		cout << "    NODE NOT NULL" << endl;

    	}

    	const uint8_t &idx = node->ieIdx;         // Reference to element

      // ---------
      // Inclusion
      // ---------

      bool isElementLeft = elements[idx]; 	// If we have the element left
      if (isElementLeft) {
      	const uint64_t inclusionSum  = pathSum  + m_S[idx]; // Inclusion sum
// 	      cout << "         Idx INC: " << (int) idx << "   IS = PS + m_S[idx]: "<< inclusionSum << " = " << pathSum << " + " << m_S[idx] << endl;

        if (pathSum == 0) {						// If firstIdx not set yet
          if (idx > firstIdxMax) {		// If the idx > firstIdxMax, then there is no way to generate
            return ub;							  // enough sets with this cardinality for a completion, so prune
          }
          firstIdx = idx+1;						// We are including an element, set it to next idx
        }
        elements.set(idx,false);      // Element no longer available

        if (cardinality == 1) {
              		cout << "  genRestCIE(" << K-1 << ", " << getBits(elements,4,elements.size()) << "...)" << endl;
              	}

        uint64_t newValue = std::max(inclusionSum,
                                     generateRestCIE(node->included, K, elements, elementsSum - m_S[idx], elementsCount-1, inclusionSum,
                                            cardinality, firstIdx, firstIdxMax, lb, partialCost, ub));
        elements.set(idx,true);           // Element available again
        ub = std::min(ub,newValue);
      } else {
//      	cout << "         Idx EXC: " << (int) idx << endl;
      }

      // ---------
      // Exclusion
      // ---------

      if (cardinality > 1 && // Don't exclude with cardinality 1. Always better to take larger numbers
					partialCost < ub &&
          pathSum < ub &&
          node->excluded != NULL) {

        uint64_t newLb;							  // For Moffitt Pruning rule
//        if (inclusionSum < ub &&      // Inclusion did not lead to a better solution
//            inclusionSum >= lb && 		// But that sum is at least as good as lower bound
//            isElementLeft) {					// And it was included (this was the bug)
//          newLb = inclusionSum + 1;	  // Then new CMin is one greater than inclusionSum
//        } else {
//          newLb = lb;								  // else use old CMin
//        }

      	newLb = lb;
//      	cout << "         Idx Exc: " << K << endl;
        uint64_t newValue =
            generateRestCIE(node->excluded, K, elements, elementsSum, elementsCount, pathSum,
                            cardinality, firstIdx, firstIdxMax, newLb, partialCost, ub);

        ub = std::min(ub,newValue);
      }
    }

    return ub;
  }




//  - Keep track of cardinality
//  - Perhaps figure out some pruning rules given the cardinality left
  uint64_t GCIEC::generateRestIE(const size_t elementIdx,      // The idx of the next element to include/exclude
                                 const int K,                  // The idx of the current bin we are filling
                                 ss::DynamicBitset &elements,  // Bitset with 1 representing elements still left
                                 const uint64_t elementsSum,   // The sum of the elements still left
                                 const size_t elementsCount,   // The count of the elements left (# of 1's in elements)
                                 const uint64_t pathSum,       // The sum of the included elements on the current path
                                 const size_t pathCardinality, // The number of included elements on the current path
                                 const size_t cardinality,     // The cardinality of sets in the tree we are searching now
                                 size_t firstIdx,              // The first idx we can include (in order to stop duplicate permutations)
                                 const size_t firstIdxMax,     // The last idx we can include so we can still make enough completions with this cardinality
                                 const uint64_t lb,            // The minimum sum of a set that could lead to a solution
                                 const uint64_t partialCost,   // The max sum used so far above this point in the tree
                                 uint64_t ub) {                // The maximum sum of a bin determined by the first set sum (k==m_k)

  	d_nodesExplored[cardinality]++;

    // ---------------------------------------------
    // All elements included or excluded (Base Case)
    // ---------------------------------------------
    if (elementIdx == boost::dynamic_bitset<>::npos || pathCardinality == cardinality) {
      if (pathSum >= lb && pathSum < ub ) {
        ub = m_partition->partition(K-1,elements,elementsSum,elementsCount,cardinality,0,std::max(partialCost,pathSum),ub);
      }
    } else {

      //if (pathSum + elementsSum >= lb) {    // If we can reach the lower bound
    	if (!canPrune(elementIdx,elements,pathSum,pathCardinality,cardinality,lb,ub)) {
        if (pathSum == 0) {                 // If firstIdx not set yet
          if (elementIdx > firstIdxMax) {   // If the idx > firstIdxMax, then there is no way to generate
            return ub;                      // enough sets with this cardinality for a completion, so prune
          }
          firstIdx = elementIdx + 1;        // We are including an element, set it to next idx
        }

        size_t nextIdx = elements.find_next(elementIdx);  // The next set bit

        //----------
        // Inclusion
        //----------
        const uint64_t newRemainingSum = elementsSum - m_S[elementIdx]; // Inclusion sum remaining
        const uint64_t inclusionSum  = pathSum  + m_S[elementIdx]; // Inclusion sum

        if (inclusionSum < ub) {
          elements[elementIdx] = false;
          uint64_t newValue = std::max(inclusionSum,
                                       generateRestIE(nextIdx, K, elements, newRemainingSum, elementsCount-1, inclusionSum,
                                      		 pathCardinality+1, cardinality, firstIdx, firstIdxMax, lb, partialCost, ub));
          elements[elementIdx] = true;
          ub = std::min(ub,newValue);
        }

        //----------
        // Exclusion
        //----------

        if (partialCost < ub &&
            pathSum < ub) {
          uint64_t newLb;               // For Moffit's pruning rule
          if (inclusionSum < ub &&      // Inclusion did not lead to a better solution
              inclusionSum >= lb) {     // But that sum is at least as good as lower bound
            newLb = inclusionSum + 1;   // Then new CMin is one greater than inclusionSum
          } else {
            newLb = lb;                 // else use old CMin
          }

          uint64_t newValue = std::max(pathSum,
                                       generateRestIE(nextIdx, K, elements, elementsSum, elementsCount, pathSum,
                                                      pathCardinality, cardinality, firstIdx, firstIdxMax, newLb, partialCost, ub));
          ub = std::min(ub,newValue);

        }
      }
    }
    return ub;
  }

  bool GCIEC::canPrune(size_t elementIdx,   // The idx of the next element to include/exclude
  		const ss::DynamicBitset &elements, 		// Bitset with 1 representing elements still left
  		const uint64_t pathSum,      					// The sum of the elements inclued so far in this tree recursion
  		const size_t pathCardinality, 				// The number of included elements on the current path
  		const size_t cardinality,     				// The cardinality of sets in the tree we are searching now
  		const uint64_t lb,       							// The minimum sum of a set that could lead to a solution
  		uint64_t ub) {                				// The maximum sum of a bin determined by the first set sum (k==m_k)

  	size_t cardinalityRemaining = cardinality-pathCardinality;	// The number of integers still to add
  	std::deque<uint64_t> queue;															// keep track of cardinalityRemianing elements
  	uint64_t sum = 0;																				// Keep running sum

  	// Look for the first cardinality remaining elements and keep the sum
  	// of the largest sum of this number of integers
  	while (elementIdx != boost::dynamic_bitset<>::npos &&		// While integers left
  				 queue.size() < cardinalityRemaining) {						// And we have not reached the min number of elements
  		sum += m_S[elementIdx];																// Add integer to sum
  		queue.push_back(m_S[elementIdx]);											// Add integer to queue
    	elementIdx = elements.find_next(elementIdx);  				// Go to next integer
  	}

  	if (pathSum + sum < lb ||										// If we can't reach lower bound with the largest integers
  			queue.size() < cardinalityRemaining) {	// Or there are not enough integers remaining to get us to cardinality
  		return true;															// We can prune
  	} else {
    	while (elementIdx != boost::dynamic_bitset<>::npos) {		// While integers left

    		sum += m_S[elementIdx];
    		sum -= queue.front();					// Add new integer and subtract largest integer
    		queue.push_back(m_S[elementIdx]);											// Remove the larget integer from the queue
    		queue.pop_front();
      	elementIdx = elements.find_next(elementIdx);  				// Go to next integer
    	}
  	}

  	return (pathSum + sum > ub);		// If the current sum + the smallest remaining elements are greater than ub, we prune

  }

  uint64_t GCIEC::generate(const int K,										// The idx of the current bin we are filling
                          ss::DynamicBitset &elements,		// Bitset with 1 representing elements still left
                          const uint64_t elementsSum,			// The sum of the elements still left
                          const size_t elementsCount,			// The count of the elements left (# of 1's in elements)
                          size_t cardinality,							// The cardinality of sets in the tree we are searching now
                          size_t firstIdx,								// The first idx we can include (in order to stop duplicate permutations)
                          uint64_t lb,										// The minimum sum of a set that could lead to a solution
                          uint64_t partialCost,						// The max sum used so far above this point in the tree
                          uint64_t ub) {									// The maximum sum of a bin determined by the first set sum (when K==m_k)


    size_t firstIdxMax;
    // Consider cardinality=4, K=9 and elementsCount = 40
    // This will succeed for cardinality=4 since 9*4 = 36 <= 40.
    // It will exit at cardinality=5 since 9*5 = 45 > 40. i.e., We can't solve
    // the problem without completing some of the bins with cardinality of 4
    while (partialCost < ub && cardinality * K <= elementsCount) {

//    	cout << "    Cardinality: " << cardinality << endl;
//    	cout << "    lb         : " << lb << endl
//    			 << "    ub         : " << ub << endl;

      if (cardinality == 1) {
      	cout << "d_timesCalled(" << getBits(elements,4,elements.size()) + ") ec: " << elementsCount << " K: " << K << endl;

      }

    	d_timesCalled[cardinality]++;

      // First do cardinality pruning
      int minRequired = ((cardinality+1) * K) - elementsCount;            // Minimum number of sets of this cardinality required

      if (minRequired > 0) {                                              // If we require any
        ss::DynamicBitset shiftedElements = elements >> firstIdx;         // Ignoring everything before firstIdx
        int shiftedElementsCount = shiftedElements.count();               // Number of 1 bits in shiftedElements
        int maxPossible = shiftedElementsCount / cardinality;             // Compute the max we can create given the first idx and the elements remaining
//        cout << "   sec : " << shiftedElementsCount << endl
//        		 << "   car : " << cardinality << endl
//        		 << "sec/car: " << maxPossible << endl;
        if (minRequired > maxPossible) {                                  // If we can't create the min required
          break;                                                          // prune
        }
//        cout << "   Min Required: " << minRequired << "  card: " << cardinality << " target: " << (cardinality*minRequired)-1 << endl;

        firstIdxMax = getFirstIdxMax(shiftedElements, firstIdx, shiftedElementsCount, (cardinality*minRequired));
      } else {
        firstIdxMax = elements.size();
      }

//      cout << "   First Idx Max: " << firstIdxMax << endl;
      // ------------------------
      // If using Cached IE Trees
      // ------------------------

      if (cardinality <= m_cachedIETrees->getMaxCardinality()) {

        const CachedIENode *node = m_cachedIETrees->getRoot(cardinality); // Choose correct root based on cardinality

        while (node != NULL && node->ieIdx < firstIdx) {
          node = node->excluded;                // Point to exclusion
        }

        ub = generateRestCIE(node, K, elements, elementsSum, elementsCount, 0,
                             cardinality, firstIdx, firstIdxMax, lb, partialCost, ub);

      } else {

        // ------------------------------------
        // If using Regular Inclusion/Exclusion
        // ------------------------------------
        size_t firstIdx = elements.find_first();

        // Create elementValues vector
        vector<uint64_t> elementValues;
        size_t elementIdx = firstIdx;
        while (elementIdx != boost::dynamic_bitset<>::npos) {

        	elementValues.push_back(m_S[firstIdx]);
        	elementIdx = elements.find_next(elementIdx);  // The next set bit
        }

        ub = generateRestIE(firstIdx, K, elements, elementsSum, elementsCount, 0,0,
                            cardinality, firstIdx, firstIdxMax, lb, partialCost, ub);
      }

      if (cardinality == 1) {
      	cout << "   K: " << K << "  GO TO NEXT CARD" << endl;
      }
      cardinality++;		// When we go to the next cardinality,


//      cout << "  pc < ub: " << partialCost << " < " << ub << " && CARD * K <= EC: " << cardinality << " * " << K << " <= " << elementsCount << endl;
      firstIdx = 0;	    // reset firstIdx
    }

    return ub;
  }


string GCIEC::cardinalityCountsToString() const {
	std::ostringstream out;
	  out << "\nCARDINALITY COUNTS" << endl;

	  out << std::setw(2)  << "i" 					<< "  "
	  		<< std::setw(12) << "CIE Subsets" 		<< "  "
	  		<< std::setw(12) << "Times Called" 		<< "  "
	  		<< std::setw(14) << "Nodes Explored" 	<< "  "
	  		<< std::setw(10) << "Ratio" << endl;

	  for (size_t i=0;i<m_S.size();i++) {

			//	  	if (d_timesCalled[i] > 0) {
			if (m_cachedIETrees->getNumSubsets(i) > 0)  {
	  		out << std::fixed << std::setprecision(2)
	  				<< std::setw(2)  << i 																	<< "  "
	  				<< std::setw(12) << m_cachedIETrees->getNumSubsets(i) 	<< "  "
	  				<< std::setw(12) << d_timesCalled[i]	 								 	<< "  "
	  				<< std::setw(14) << d_nodesExplored[i] 								<< "  "
	  				<< std::setw(10) << (double) d_nodesExplored[i] / (double) d_timesCalled[i] << endl;
	  	}
	  }
	  return out.str();
}

}


