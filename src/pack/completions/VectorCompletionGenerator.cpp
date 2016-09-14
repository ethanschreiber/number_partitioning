/*
 * VectorCompletion.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#include "VectorCompletionGenerator.hpp"

namespace ss {
bool VectorCompletionGenerator::next(SetNodeVector &node) {
  if (m_idx < m_completions.size()) {
    node = m_completions[m_idx];
    m_idx++;
    return true;
  } else {
    return false;
  }
}


} /* namespace ss */
