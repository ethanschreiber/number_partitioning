/*
 * Completion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef COMPLETION_GENERATOR_HPP_
#define COMPLETION_GENERATOR_HPP_
#include "../PackingUtils.hpp"
#include "../../Heap.hpp"

namespace ss {
class CompletionGenerator {
public :
  virtual bool next(SetNodeVector &node) = 0;
  virtual ~CompletionGenerator() {};
};

} // end namespace

#endif /* COMPLETION_GENERATOR_HPP_ */
