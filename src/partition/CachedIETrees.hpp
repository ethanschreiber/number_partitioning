/*
 * MoffittCachedIE.hpp
 *
 *  Created on: Dec 24, 2013
 *      Author: ethan
 */

#ifndef CACHEDIE_HPP_
#define CACHEDIE_HPP_

#include "../ss/SubsetSum.hpp"
#include "CachedSetsBuffered.hpp"

namespace partition {

// ============================================================================
//
// CachedIENode - Struct used by CachedIETrees classes.
//
// ============================================================================
struct CachedIENode {

  CachedIENode *included;
  CachedIENode *excluded;
  uint8_t ieIdx;  // The index of the element to include/exclude


  CachedIENode() : included(NULL), excluded(NULL), ieIdx(0) {
  }

  CachedIENode(const uint8_t ieIdx) : included(NULL), excluded(NULL), ieIdx(ieIdx) {
  }

  CachedIENode(const uint8_t ieIdx, CachedIENode *included, CachedIENode *excluded) :
    included(included), excluded(excluded),ieIdx(ieIdx) {

  }
};

// ============================================================================
//
// C a c h e d I E T r e e s
//
// CachedIETrees is abstract, getCachedSetsBuffered() is pure virtual
//
// ============================================================================

class CachedIETrees {

protected :
  const vector<uint64_t> m_S; // The input elements
  uint64_t m_elementsSum;     // The sum of all elements in S
  int m_K;                    // The total number of bins in the problem
  CachedIENode  **m_roots;    // dynamic array of root pointers to prefix trees

  vector<size_t> d_numSubsetsPerCard; // Debugging, keep track of number of subsets of each cardinality

  // Add node to proper cached IE tree
  void push_back(const ss::SetNodeBitset &node);

  // Clear an IE tree below a node (i.e., for one cardinality tree)
  void clear(CachedIENode *node);

  // Clear all Cached IE trees
  void clear();

  virtual CachedSetsBuffered* getCachedSetsBuffered() const = 0;
public :
  CachedIETrees(const vector<uint64_t> &S, int K);
  virtual ~CachedIETrees();

  // Add sets to the CIE Tree
  // Return true if succesful, false if none left to add
  bool addSetsToCache();

  const CachedIENode* getRoot(size_t cardinality) const;

  // Debugging, report number of subsets of a particular cardinality
  size_t getNumSubsets(size_t cardinality);
  // Return the minimal cardinality of any set
  // stored in this data structure.
  // returns N if no sets are stored
  size_t getMinCardinality() const;

  // Return the maximal cardinality of any set
  // stored in this data structure.
  // returns N if no sets are stored
  size_t getMaxCardinality() const;

  // Get the max subset in range
    SetVectorIt getMax() const {
      return getCachedSetsBuffered()->getMax();
    }

    // Get the subsets with sum below perfect in range
    SetIteratorPair getMinRange() const {
      return getCachedSetsBuffered()->getMinRange();
    }

  // ===================
  // Debugging Functions
  // ===================
    // Number of times Schroeppel and Shamir was called
    size_t getSsCalls() const {
      return getCachedSetsBuffered()->getGenerationCalls();
    }

    // Amount of time spent running Schroeppel and Shamir
    double getSsTime() const {
      return getCachedSetsBuffered()->getGenerationTime();
    }

#if PAPER_DEBUG == 1
protected :
  void toDot(std::ostringstream &out,  string path,CachedIENode  *node);
public :
  // Debugging functions
  string toDot();
  string printSets(uint64_t min, uint64_t max);
#endif

};

// ============================================================================
//
// S u b c l a s s e s
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// C a c h e d I E T r e e s A l l C a r d i n a l i t y S S
//
// The cached Iterative Weakening trees for CIW
// This uses Schroeppel and Shamir to generate "numSets" subsets with sum
// greater than or equal to perfect at a time as well as the corresponding subsets
// with sum below perfect
//
// ----------------------------------------------------------------------------
class CachedIETreesAllCardinalitySS : public CachedIETrees {
protected :
  CachedSetsBufferedSS* m_cs;
  CachedSetsBuffered* getCachedSetsBuffered() const {
    return m_cs;
  }
public:
  CachedIETreesAllCardinalitySS(const vector<uint64_t> &S, int K, uint64_t upperBound, size_t numSets);
	virtual ~CachedIETreesAllCardinalitySS();
};

// ----------------------------------------------------------------------------
//
// C a c h e d I E T r e e s L o w C a r d i n a l i t y H S
//
// ----------------------------------------------------------------------------
class CachedIETreesLowCardinalityHS: public CachedIETrees {
protected :
  CachedSetsBufferedLowCardinalityHS* m_cs;
  CachedSetsBuffered* getCachedSetsBuffered() const {
    return m_cs;
  }
public:

  CachedIETreesLowCardinalityHS(const vector<uint64_t> &S, const int K, const uint64_t upperBound,
                                const size_t numSets, const size_t maxCardinality);

	virtual ~CachedIETreesLowCardinalityHS();
};


} /* namespace partition */

#endif /* CACHEDIE_HPP_ */
