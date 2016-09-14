/*
 * VectorCompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef OPSCOMPLETION_GENERATOR_HPP_
#define OPSCOMPLETION_GENERATOR_HPP_

#include "Completion_Generator.hpp"
#include "../../ss/OrderedPowerSet.hpp"

#include <vector>
using std::vector;

namespace ss {

class OPSCompletionGenerator: public ss::CompletionGenerator {
private :
  const int m_N;
  const uint64_t *m_S;
  const uint64_t m_lower;
  const uint64_t m_upper;
  PersisistentIndexDeque<SetNodeBitset> m_L;
  OPSMinHeap m_nextHeap;
public :
  OPSCompletionGenerator(const uint64_t *S, const int N,
      const uint64_t lower, const uint64_t upper, deque<SetNodeBitset> &LParam);
  bool next(SetNodeVector &node);

  ~OPSCompletionGenerator();
};

} /* namespace ss */
#endif /* OPSCOMPLETION_GENERATOR_HPP_ */
