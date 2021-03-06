#include "../../PackingUtils.hpp"
#include <stdint.h>
#include <iostream>
#include <vector>
#include <deque>
#include <stack>

using namespace std;

const int SIZE=3;
const int64_t m_elements[SIZE] = {3,2,1};
SetNodeVector m_included;
SetNodeVector m_excluded;

struct IENode {
  const int64_t *unassigned;    // Pointer to first unassigned element in array
  const int64_t sumUnassigned;  // Sum of unassigned elements
  const int64_t subsetSum;      // Subset sum of included elements
  const int64_t lBound;         // The lower bound
  const int     numUnassigned;  // Number of unassigned elements
  const int     numIncluded;    // The number included before taking action on this node
  const int     numExcluded;    // The number excluded before taking action on this node
  const bool    isInclusion;    // True if inclusion node, false if exclusion.

  IENode(const int64_t *unassigned, const int64_t sumUnassigned, const int64_t subsetSum,
         const int numUnassigned, const int numIncluded, const int numExcluded, const int64_t lBound, const bool isInclusion)
  : unassigned(unassigned), sumUnassigned(sumUnassigned), subsetSum(subsetSum),
    numUnassigned(numUnassigned), numIncluded(numIncluded), numExcluded(numExcluded),
    lBound(lBound),isInclusion(isInclusion){}
};

stack<IENode> m_stack;


int main() {

  int64_t sumElements = 0;
  for (int i=0;i<SIZE;i++) {
    sumElements+= m_elements[i];
  }
  int64_t sumUnassigned = sumElements - m_elements[0];
  int numUnassigned     = SIZE - 1;

  m_stack.push(IENode(m_elements, sumUnassigned, 0            , numUnassigned,0,1, -1, false));
  m_stack.push(IENode(m_elements, sumUnassigned, m_elements[0], numUnassigned,1,0, -1, true));

  while (!m_stack.empty()) {
    const int64_t* unassigned    = m_stack.top().unassigned;    // Copy unassigned pointer.
    int64_t        sumUnassigned = m_stack.top().sumUnassigned; // The sum of the unassigned elements.
    int            numUnassigned = m_stack.top().numUnassigned; // The number of unassigned elements.
    int            numIncluded   = m_stack.top().numIncluded;   // The number included.
    int            numExcluded   = m_stack.top().numExcluded;   // The number excluded.
    const int64_t  subsetSum     = m_stack.top().subsetSum;     // The subset sum of the included elements.
    const int64_t  lBound        = m_stack.top().lBound;        // The lower bound.
    const bool     isInclusion   = m_stack.top().isInclusion;   // Is this an inclusion node?

    m_stack.pop();      // Pop node

    int            firstIdx      = 0;                           // Index of largest element that fits in bin.

    // --------------------------------
    // Modify m_included and m_excluded
    // --------------------------------

    m_included.pop_back(m_included.cardinality() - numIncluded);
    m_excluded.pop_back(m_excluded.cardinality() - numExcluded);
    if (isInclusion) {
      m_included.push_back(unassigned[0]);
    } else {
      if (numExcluded && m_excluded.cardinality() == numExcluded) {
        m_excluded.pop_back();
      }
      m_excluded.push_back(unassigned[0]);
    }

    // -----------------------------------------
    // Check for leaf nodes and modify the stack
    // -----------------------------------------
    if (numUnassigned == 0) {
      cout << "Included: " << m_included.toString() << " " << subsetSum << endl;
      cout << "Excluded: " << m_excluded.toString() << endl << endl;
    } else {

      m_stack.push(IENode(unassigned+1, sumUnassigned-unassigned[1], subsetSum                  , numUnassigned-1, numIncluded  , numExcluded+1, -1, false));

      // Never allow an exclusion followed by an inclusion of a duplicate value
      if (m_excluded.empty() || unassigned[1] != m_excluded.back()) {
        m_stack.push(IENode(unassigned+1, sumUnassigned              , subsetSum + unassigned[1], numUnassigned-1, numIncluded+1, numExcluded  , -1, true));
      }
    }
  }
}
