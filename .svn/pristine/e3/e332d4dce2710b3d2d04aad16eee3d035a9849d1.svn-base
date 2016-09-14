/*
 * Moffitt.cpp
 *
 *  Created on: Jul 31, 2013
 *      Author: ethan
 */


#include "GenerateIEBitset.hpp"
#include "PartitionUtils.hpp"

namespace partition {



GenerateIEBitset::GenerateIEBitset(Partition* partition, const vector<uint64_t>& S, int K)
: m_partition(partition), m_S(S.begin(),S.end()), m_K(K)
{
}


uint64_t GenerateIEBitset::generateRest(const int K, ss::DynamicBitset &elements, const uint64_t elementsSum,
                      									const size_t elementsCount, const uint64_t includedSum,
                      									const size_t elementIdx, size_t cardinality,
                      									const uint64_t lb, uint64_t partialCost, uint64_t ub) {



	// ---------------------------------------------
  // All elements included or excluded (Base Case)
  // ---------------------------------------------
  if (elementIdx == boost::dynamic_bitset<>::npos) {

  	if (includedSum >= lb && includedSum < ub) {

    	ub = m_partition->partition(K-1,elements,elementsSum,elementsCount,cardinality,0,std::max(partialCost,includedSum),ub);
    }
  } else {

    if (includedSum + elementsSum >= lb) {   // If we can reach the lower bound

      size_t nextIdx = elements.find_next(elementIdx);  // The next set bit

      //----------
      // Inclusion
      //----------
      const uint64_t newRemainingSum = elementsSum - m_S[elementIdx]; // Inclusion sum remaining
      const uint64_t newIncludedSum  = includedSum  + m_S[elementIdx]; // Inclusion sum

      if (newIncludedSum < ub) {
        elements[elementIdx] = false;
        uint64_t newValue = std::max(newIncludedSum,
                                     generateRest(K,elements,newRemainingSum, elementsCount-1, newIncludedSum,
                                    		 	 	 	 	 	nextIdx, cardinality, lb, partialCost, ub));
        elements[elementIdx] = true;
        ub = std::min(ub,newValue);
      }

      //----------
      // Exclusion
      //----------

      if (partialCost < ub &&
          includedSum < ub) {
        uint64_t newLb;               // For Moffit's pruning rule
        if (newIncludedSum < ub && newIncludedSum >= lb) {
          newLb = newIncludedSum + 1;
        } else {
          newLb = lb;
        }

        uint64_t newValue = std::max(includedSum,
                                     generateRest(K, elements, elementsSum, elementsCount, includedSum,
                                    		 	 	 	 	 	nextIdx, cardinality, newLb, partialCost, ub));
        ub = std::min(ub,newValue);

      }
    }
  }
  return ub;
}

uint64_t GenerateIEBitset::generate(const int K,
                  								  boost::dynamic_bitset<> &elements,
                  								  const uint64_t elementsSum, const size_t elementsCount,
                  								  size_t cardinality, const uint64_t lb, uint64_t partialCost, uint64_t ub) {


  if (cardinality >= 6) {
    cout << "K: " << std::setw(2) << K << " " << "Card: " << cardinality << "  EC: " << elementsCount << endl;
  }

	size_t firstIdx = elements.find_first();
  size_t nextIdx  = elements.find_next(firstIdx);

  elements[firstIdx] = false;
  uint64_t val = generateRest(K, elements, elementsSum-m_S[firstIdx], elementsCount-1, m_S[firstIdx],
  														nextIdx, cardinality, lb, partialCost, ub);

  elements[firstIdx] = true;
  return val;

}


} // End Namespace
