/*
 * Partition.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: ethan
 */

#include "Partition.hpp"


namespace partition {
  Partition::Partition(const vector<uint64_t>& S, const int K)
  : m_S(S.begin(),S.end()), m_K(K), m_firstCount(0) {

    d_binCounts.resize(K + 1);

    for (int i=0;i<K + 1;i++) {
      d_binCounts[i] = 0;
    }
  }

  Partition::~Partition() {

  }

  size_t Partition::getFirstCount() const {
    return m_firstCount;
  }

}
