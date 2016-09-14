/*
 * SubsetSum.hpp
 *
 *  Created on: Aug 29, 2012
 *      Author: ethan
 */

#ifndef SUBSETSUM_HPP_
#define SUBSETSUM_HPP_

#include "../pack/PackingUtils.hpp"

#include <vector>
#include <deque>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
using std::vector;
using std::deque;

namespace ss {

typedef boost::dynamic_bitset<uint64_t> DynamicBitset;

// ============================================================================
// SetNodeBitset Class - This is the SetNode class that stores elements as a
// bitset. Not to be confused with SetNodeVector.
// ============================================================================
struct SetNodeBitset {
private :
  uint64_t m_sum;
  DynamicBitset m_set;
  uint8_t m_lastSetIdx;
  uint8_t m_cardinality;


public :
  static const DynamicBitset::size_type npos = DynamicBitset::npos;
  SetNodeBitset() {}
  SetNodeBitset(uint64_t sum0,const DynamicBitset &set0) : m_sum(sum0), m_set(set0), m_cardinality(set0.count()) {}
  SetNodeBitset(size_t size) : m_sum(0), m_set(size), m_cardinality(0) {}

  // This constructor only works if the bits set in set 1
  // are disjoint from the bits set in set 2
  SetNodeBitset(const SetNodeBitset &set1,const SetNodeBitset &set2) :
  m_sum(set1.sum() + set2.sum()),
  m_set(set1.set() | set2.set()),
  m_cardinality(set1.cardinality() + set2.cardinality()){ }
  void set(size_t index, uint64_t value) {
    m_set[index] = true;
    m_sum += value;
    m_cardinality++;
    m_lastSetIdx = index;
  }

  void unset(size_t index, uint64_t value) {
    m_set[index] = false;
    m_sum -= value;
    m_cardinality--;
  }


  uint64_t sum() const {
    return m_sum;
  }

  uint8_t cardinality() const {
    return m_cardinality;
  }


  uint64_t reverseSum(uint64_t maxValue) const {
    return maxValue - m_sum;
  }

  const DynamicBitset& set() const {
    return m_set;
  }

  string toString(const uint64_t *S) const {
    std::ostringstream out;
    out << "[card = " << m_set.count() << " members = {";
    size_t idx = m_set.find_first();
    bool isFirst = true;
    while (idx != npos){
      if (isFirst) {
        isFirst = false;
      } else {
        out << ",";
      }

      out << S[idx];
      idx = m_set.find_next(idx);
    }
    out << "} sum = " << m_sum << "]";
    return out.str();
  }

  string toString() const {
    std::ostringstream out;
    out << "[sum = " << sum() << " card = " << std::setfill('0') << std::setw(2) << m_set.count() << " members = {";
    size_t idx = m_set.find_first();
    bool isFirst = true;
    while (idx != npos){
      if (isFirst) {
        isFirst = false;
      } else {
        out << ",";
      }

      out << std::setfill('0')<< std::setw(2) << idx;
      idx = m_set.find_next(idx);
    }
    out << "} ]";
    return out.str();
  }


  DynamicBitset::size_type max()  const {
  	cout << "ETHAN, look at max. This should only be used in OrderedPowerSet.cpp and OPSCompletionGenerator.cpp\n"
  			<< "This doesn't seem to be correct as m_lastSetIdx never changes\n\n";
//    DynamicBitset::size_type idx = m_set.find_first();
//    DynamicBitset::size_type lastIdx;
//
//    do {
//      lastIdx = idx;
//      idx = m_set.find_next(idx);
//    } while (idx != npos);
//    if (lastIdx != m_lastSetIdx) {
//    	std::cerr << endl << lastIdx << " != " << m_lastSetIdx << endl;
//    	exit(0);
//    }
//    return lastIdx;
  	return m_lastSetIdx;
  }

  bool operator[](size_t idx) const{
    return m_set[idx];
  }


  friend inline bool operator== (const SetNodeBitset &node1, const SetNodeBitset &node2);
  friend inline bool operator< (const SetNodeBitset &node1, const SetNodeBitset &node2);
  friend inline bool operator> (const SetNodeBitset &node1, const SetNodeBitset &node2);

};


inline bool operator== (const SetNodeBitset &node1, const SetNodeBitset &node2) {
    return (node1.m_sum == node2.m_sum) && (node1.m_set == node2.m_set);
}

inline bool operator< (const SetNodeBitset &node1, const SetNodeBitset &node2) {
    return node1.m_sum < node2.m_sum;
}

inline bool operator> (const SetNodeBitset &node1, const SetNodeBitset &node2) {
    return node1.m_sum > node2.m_sum;
}

class SetComparatorAscending
{
public:
    bool operator()(const SetNodeBitset& node1, const SetNodeBitset& node2)
    {
        return node1 > node2;
    }
};

// -------
// General
// -------

const uint8_t UNSET_MAX_CARD=std::numeric_limits<uint8_t>::max();
void generateAllSets(const uint64_t S[], const int N,
             vector<SetNodeBitset> &sets,
             int first, int last,
             bool sortAscending = true, uint8_t maxCardinality=UNSET_MAX_CARD);

// Helper function, provides missing variables
// also sorts sets before returning
void generateAllSums(const uint64_t S[],
             	 	 	 	 vector<uint64_t> &sums,
             	 	 	 	 int first, int last,
             	 	 	 	 bool sortAscending = true);

// --------------------------
// Schroeppel and Shamir Half
// --------------------------
size_t generateSumsSSHalf(const uint64_t S[], const int n, const uint64_t lower, const uint64_t upper,
             vector<uint64_t> &allsums);


size_t generateSetsSSHalf(const uint64_t S[], const int n,
                          const uint64_t lower, const uint64_t upper,
                          vector<SetNodeBitset> &allSets);

// -----------------
// Ordered Power Set
// -----------------

size_t generateSetsOPS(const uint64_t S[], const int n,
                       const uint64_t lower, const uint64_t upper,
                       deque<SetNodeBitset> &allSets);

size_t generateSetsOPS(const uint64_t S[], const int n,
                       const uint64_t lower, const uint64_t upper,
                       deque<SetNodeVector> &allSets);


size_t generateSetsOPSReverse(const uint64_t S[], const int N,
                              const uint64_t lower, const uint64_t upper,
                              deque<SetNodeBitset> &LParam);

} // end namespace
#endif /* SSALLSETS_HPP_ */
