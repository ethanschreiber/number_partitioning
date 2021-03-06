/*
 * BufferedCompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef BUFFEREDCOMPLETION_GENERATOR_HPP_
#define BUFFEREDCOMPLETION_GENERATOR_HPP_

#include "Completion_Generator.hpp"

#include <vector>
using std::vector;

namespace ss {

class BufferedCompletionGenerator: public ss::CompletionGenerator {
private :
  CompletionGenerator *m_generator;
  size_t m_idx;
  vector<SetNodeVector> m_buffer;
  size_t m_bufferSize;
  bool m_isDone;
  bool m_isClassicSort;		// Should we use the classic Korf sort?
protected :
  bool fillBuffer() {
    m_idx = 0;
    m_buffer.clear(); // Clear the buffer

    if (!m_isDone) {
      // Fill the buffer
      SetNodeVector node;
      while (m_buffer.size() < m_bufferSize && m_generator->next(node)) {

        m_buffer.push_back(node);
      }
      m_isDone = (m_buffer.size() < m_bufferSize);

      if (m_isClassicSort) {
    	  std::sort(m_buffer.begin(),m_buffer.end(),SetNodeCardinalityComparatorClassic());
      } else {
    	  std::sort(m_buffer.begin(),m_buffer.end(),SetNodeCardinalityComparator());
      }
    }
    return !m_buffer.empty();
  }
public :
  BufferedCompletionGenerator(CompletionGenerator *generator, size_t bufferSize,bool isClassicSort)
    : m_generator(generator),m_idx(0), m_bufferSize(bufferSize), m_isDone(false), m_isClassicSort(isClassicSort) {

    m_buffer.reserve(bufferSize);  // Reserve memory for buffer

  }

  bool next(SetNodeVector &node);

  ~BufferedCompletionGenerator() {
    delete m_generator;
  }
};

} /* namespace ss */
#endif /* BUFFEREDCOMPLETION_GENERATOR_HPP_ */
