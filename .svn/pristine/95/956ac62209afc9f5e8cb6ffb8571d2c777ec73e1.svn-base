/*
 * IECompletion.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: ethan
 */

#ifndef IECOMPLETION_GENERATOR_HPP_
#define IECOMPLETION_GENERATOR_HPP_

#include "Completion_Generator.hpp"

#include <vector>
#include <deque>
#include <stack>
using std::vector;
using std::deque;
using std::stack;

namespace ss {
struct IENode {
  const int         assignedIdx;    // The idx of the element to be included or excluded (unassigned starts at assignedIdx+1)
  const uint64_t    sumUnassigned;  // Sum of unassigned elements (starting at assignedIdx+1)
  const uint64_t    subsetSum;      // Subset sum of included elements (including assignedIdx if included)
  const int         numIncluded;    // The number included before taking action on this node
  const int         numExcluded;    // The number excluded before taking action on this node
  const bool        isInclusion;    // True if inclusion node, false if exclusion.

  const uint64_t lowerBound;     // The lower bound for this node. Because of moffitt, this can be higher
                                // than the global lower bound
  IENode(const int assignedIdx, const uint64_t sumUnassigned, const uint64_t subsetSum,
         const int numIncluded, const int numExcluded, const bool isInclusion, const uint64_t lowerBound)
  : assignedIdx(assignedIdx), sumUnassigned(sumUnassigned), subsetSum(subsetSum),
    numIncluded(numIncluded), numExcluded(numExcluded)    , isInclusion(isInclusion), lowerBound(lowerBound) {
  }
};


class IECompletionGenerator: public ss::CompletionGenerator {
private :


  const uint64_t m_lBound;   // The lower bound for completions
  const uint64_t m_uBound;   // The upper bound for completions

  const int m_numElements;  // The number of elements in the m_elements array
  uint64_t m_elements[MAXN]; // The array of all elements we can use to create completions

  SetNodeVector m_included;       // While searching, contains the included elements at any node in tree.
  SetNodeVector m_excluded;       // While searching, contains the excluded elements at any node in tree.

  stack<IENode> m_stack;    // The stack of nodes in support of the DFS branch and bound search

  bool m_checkDominance;    // Should we check for dominance?
public :
  static size_t s_tmpCount;
  IECompletionGenerator(const uint64_t lBound,          // the lower bound for the subset sum
                        const uint64_t uBound,          // the upper bound for the subset sum
                        const uint64_t elements[MAXN],// array of unassigned elements
                        const int numElements,       // the number of unassigned elements
                        const uint64_t sumElements,     // the sum of the elements
                        bool checkDominance)

: m_lBound(lBound), m_uBound(uBound), m_numElements(numElements), m_checkDominance(checkDominance)
{
    memcpy(m_elements,elements,sizeof(uint64_t) * numElements);

    // Put first two nodes on stack
    uint64_t sumUnassigned = sumElements;

    int newAssignedIdx = 0;                      // index of largest element that fits in bin
    while (newAssignedIdx < numElements &&       // unssigned elements left and
           elements[newAssignedIdx] > uBound) {  // element larger than capacity
      sumUnassigned -= elements[newAssignedIdx]; // subtract it from sum of remaining elements
      newAssignedIdx++;                          // Start at next element
    }

    const uint64_t newAssigned   = elements[newAssignedIdx];
    const uint64_t newSumUnassigned = sumUnassigned - newAssigned;
//
//    cout << "-------------------------" << endl
//         << "Init LB      : " << lBound << endl
//         << "Init UB      : " << uBound << endl
//         << "Assigned     : " << newAssigned << endl
//         << "Assigned idx : " << newAssignedIdx << endl
//         << "Sum Unassigne: " << newSumUnassigned << endl
//         << "-------------------------" << endl << endl;

    m_stack.push(IENode(newAssignedIdx, newSumUnassigned, 0          , 0, 1, false,m_lBound));
    m_stack.push(IENode(newAssignedIdx, newSumUnassigned, newAssigned, 1, 0, true ,m_lBound));

}

  void setExcluded(SetNodeVector & excluded) const {

    size_t includedIdx = 0;
    int elementIdx;
    const vector<uint64_t> & included = m_included.elements();

    // Put all elements not in included in the excluded set
    for (elementIdx=0;elementIdx<m_numElements && includedIdx < included.size();elementIdx++) {
      if (m_elements[elementIdx] != included.at(includedIdx)) {
        excluded.push_back(m_elements[elementIdx]);
      } else {
        includedIdx++;
      }
    }

    // Put remaining elements in
    while (elementIdx < m_numElements) {
      excluded.push_back(m_elements[elementIdx]);
      elementIdx++;
    }
  }

  bool next(SetNodeVector &node);

  ~IECompletionGenerator() {}

};

} /* namespace ss */
#endif /* IECOMPLETIONGENERATOR_HPP_ */
