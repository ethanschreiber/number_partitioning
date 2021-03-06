/*
 * CachedSets.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: ethan
 */

#include "CachedSetsBuffered.hpp"
#include <numeric> // for accumulate

namespace partition {


	CachedSetsBuffered::CachedSetsBuffered(const uint64_t S[], const int N, const int K, uint64_t upperBound, size_t numSets) :
			m_elementsSum(std::accumulate(S,S+N,(uint64_t) 0)), m_K(K),
			m_smallSetIdx(0), m_largeSetIdx(0), m_generationTime(0), m_generationCalls(0)
	{
		//generateSetsSS();
	}

	CachedSetsBuffered::~CachedSetsBuffered() {
	}

	void CachedSetsBuffered::generateSets() {
		SimpleTimer timer;

		generateSets(m_smallSets,m_largeSets);    // Fill in m_smallSets and m_largeSets
		m_generationTime += timer.timeElapsed();					// Keep track of time used
		m_generationCalls++;

		m_smallSetIdx = 0;    // Reset indices
		m_largeSetIdx = 0;
	}

	bool CachedSetsBuffered::next() {
		if (m_largeSetIdx == m_largeSets.size()) {      // If we have run out of large sets
			if (!m_smallSetIdx == m_smallSets.size()) {   // We better have run out of small sets
				cout << "ERROR: Small Sets Not Empty: "
						 << m_smallSetIdx << " != " << m_smallSets.size() << endl;
				exit(0);
			}

			generateSets();

			if (m_largeSets.size() == 0) {  // If no more sets
			  return false;
			}
		}

		m_currentMax = m_largeSets.begin() + m_largeSetIdx;
		m_largeSetIdx++;

		// Compute the minimum value that can be included in a solution

		uint64_t minValue = computeCMin(m_elementsSum,m_K,m_currentMax->sum()+1);

		m_currentMinBegin = m_smallSets.begin() + m_smallSetIdx;
		while (m_smallSetIdx < m_smallSets.size() &&    // While sets left
					 m_smallSets[m_smallSetIdx].sum() >= minValue) {   // And sum is > than the minValue
			m_smallSetIdx++;
		}

		m_currentMinEnd = m_smallSets.begin() + m_smallSetIdx;
		return true;
	}

  // ---------------
  // Debugging Calls
  // ---------------

  // Number of times Schroeppel and Shamir was called
  size_t CachedSetsBuffered::getGenerationCalls() const {
    return m_generationCalls;
  }

  // Amount of time spent running Schroeppel and Shamir
  double CachedSetsBuffered::getGenerationTime() const {
    return m_generationTime;
  }


	// ============================================================================
	//
	// S u b c l a s s e s
	//
	// ============================================================================

	// ----------------------------------------------------------------------------
	//
	// C a c h e d S e t s B u f f e r e d S S
	//
	// ----------------------------------------------------------------------------

	CachedSetsBufferedSS::CachedSetsBufferedSS(const uint64_t S[], const int N, const int K, uint64_t upperBound, size_t numSets) :
			CachedSetsBuffered(S,N,K,upperBound,numSets), m_ssms(S,N,K,numSets,upperBound)
	{
	}

	void CachedSetsBufferedSS::generateSets(SetVector& smallSets, SetVector& largeSets) {
		m_ssms.generateSetsSS(m_smallSets,m_largeSets); // Fill in m_smallSets and m_largeSets
	}

	// ----------------------------------------------------------------------------
	//
	// C a c h e d S e t s B u f f e r e d L o w C a r d i n a l i t y H S
	//
	// ----------------------------------------------------------------------------

	CachedSetsBufferedLowCardinalityHS::CachedSetsBufferedLowCardinalityHS(const uint64_t S[], const int N, const int K, uint64_t upperBound,
					size_t numSets, const size_t maxCardinality) :
    CachedSetsBuffered(S,N,K,upperBound,numSets),
    //m_ssms(S,N,K,numSets,upperBound,maxCardinality)
    m_hsms(S,N,K,numSets,upperBound,maxCardinality)
	{
	}

	void CachedSetsBufferedLowCardinalityHS::generateSets(SetVector& smallSets, SetVector& largeSets) {
//	  m_ssms.generateSetsSS(m_smallSets,m_largeSets);    // Fill in m_smallSets and m_largeSets
	  m_hsms.generateSetsHS(m_smallSets,m_largeSets);    // Fill in m_smallSets and m_largeSets
	}
}
