/*
 * PartitionAlgorithm.h
 *
 *  Created on: Jul 3, 2013
 *      Author: ethan
 */

#ifndef PARTITIONALGORITHM_H_
#define PARTITIONALGORITHM_H_
#include <stdint.h>
#include "../../PackingUtils.hpp"
namespace ss {

class PartitionAlgorithm {
private :
protected :
  const uint64_t *S;
  const int N;
  const uint64_t SUM;

  uint64_t m_nodeCount;
  uint64_t m_best;

public :
  PartitionAlgorithm(const PackingProblem &p);
  virtual uint64_t execute() = 0;
};

} /* namespace ss */
#endif /* PARTITIONALGORITHM_H_ */
