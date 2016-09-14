/*
 * VectorCompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef VECTORCOMPLETION_GENERATOR_HPP_
#define VECTORCOMPLETION_GENERATOR_HPP_

#include "Completion_Generator.hpp"

#include <vector>
using std::vector;

namespace ss {

class VectorCompletionGenerator: public ss::CompletionGenerator {
private :
  const vector<SetNodeVector> &m_completions;
  size_t m_idx;
public :
  VectorCompletionGenerator(const vector<SetNodeVector> &completions)
    : m_completions(completions), m_idx(0) {}

  bool next(SetNodeVector &node);

  ~VectorCompletionGenerator() {}
};

} /* namespace ss */
#endif /* VECTORCOMPLETION_GENERATOR_HPP_ */
