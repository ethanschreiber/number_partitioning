/*
 * VectorCompletion.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#include "BufferedCompletionGenerator.hpp"

namespace ss {
bool BufferedCompletionGenerator::next(SetNodeVector &node) {

  if (m_idx < m_buffer.size()) {
    node = m_buffer[m_idx];
    m_idx++;
    return true;
  } else if (fillBuffer()){ // This resets m_idx and refills m_buffer
    node = m_buffer[m_idx];
    m_idx++;
    return true;
  } else {
    return false;
  }
}


} /* namespace ss */
