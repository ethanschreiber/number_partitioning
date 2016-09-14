/*
 * PartitionOPS.hpp
 *
 *  Created on: Jul 3, 2013
 *      Author: ethan
 */

#ifndef PARTITIONOPS_HPP_
#define PARTITIONOPS_HPP_

#include "PartitionAlgorithm.hpp"
#include "../OrderedPowerSet.hpp"

namespace ss {

struct OPSStruct {
  uint64_t sum;
  uint8_t maxIdx;
};

class PartitionOPS: public ss::PartitionAlgorithm {
private :
  PersisistentIndexDeque<OPSStruct> m_L;
  deque<OPSStruct> &m_d;

  OPSHeapNode
public:
  PartitionOPS::PartitionOPS(const PackingProblem &p);
  virtual ~PartitionOPS();

  uint64_t execute() = 0;

};

} /* namespace ss */
#endif /* PARTITIONOPS_HPP_ */
