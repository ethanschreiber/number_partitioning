/*
 * IECompletion.cpp
 *
 *  Created on: Dec 31, 2012
 *      Author: ethan
 */

#include "IECompletionGenerator.hpp"
#include "../BinCompletionUtils.hpp"
namespace ss {

bool IECompletionGenerator::next(SetNodeVector &node) {

//  The problem has to do with when to include/exclude elements. Check calculation of firstIdx
//  and when m_included/m_excluded are modified.

  // Main Loop

  while (!m_stack.empty()) {  // While there are more nodes to search

    // ===========================================================
    // Copy the elements from the top of the stack then pop node.
    // ===========================================================
    const int       assignedIdx   = m_stack.top().assignedIdx;   // Copy unassigned pointer.
    uint64_t        sumUnassigned = m_stack.top().sumUnassigned; // The sum of the unassigned elements.
    int             numIncluded   = m_stack.top().numIncluded;   // The number included.
    int             numExcluded   = m_stack.top().numExcluded;   // The number excluded.
    const uint64_t  subsetSum     = m_stack.top().subsetSum;     // The subset sum of the included elements.
    const bool      isInclusion   = m_stack.top().isInclusion;   // Is this an inclusion node?
    const uint64_t  lowerBound    = m_stack.top().lowerBound;    // The lower bound
    m_stack.pop();       // Pop node

    // ================================================================
    // Modify m_included and m_excluded for current node
    // ================================================================

    // Remove elements assigned to be included or excluded that are
    // no longer assigned for this node. Assign the new node
    {
      const uint64_t &assigned     = m_elements[assignedIdx];     // The element to be included/excluded
      m_included.pop_back(m_included.cardinality() - numIncluded);
      m_excluded.pop_back(m_excluded.cardinality() - numExcluded);
      if (isInclusion) {
        m_included.push_back(assigned);             // Include assigned element
      } else {
        if (numExcluded && m_excluded.cardinality() == numExcluded) {
          m_excluded.pop_back();
        }

        m_excluded.push_back(assigned);             // Exclude assigned element
      }
    }

//    cout << "Included        : " << m_included.toString() << endl
//         << "Excluded        : " << m_excluded.toString() << endl
//         << "Unassigned      : {" << join(m_elements+assignedIdx+1,m_numElements-(assignedIdx+1),",") << "}" << endl << endl;

    // ================================================================
    // Find the largest element that fits with the current subset
    // ================================================================
    int newAssignedIdx = assignedIdx + 1;
    while (newAssignedIdx < m_numElements &&                    // Elements left &&
           subsetSum + m_elements[newAssignedIdx] > m_uBound) { // largest can't be added to current subset
      sumUnassigned -= m_elements[newAssignedIdx];              // subtract it from sum of remaining elements
      newAssignedIdx++;                                         // skip this element
    }
    const uint64_t &currentElement = m_elements[newAssignedIdx]; // value of next largest element


    // ================================================================
    // Check Termination Conditions
    // ================================================================

    // No Elements left, add if the elements in so far are greater than lower bound
    if (newAssignedIdx == m_numElements) {                       // no elements left
      if (subsetSum >= lowerBound) {              // check to see if lower bound achieved

        if (!m_checkDominance ||
            !bp::isDominated(m_included.elements(),m_included.getSum(),
                             m_excluded.elements(),m_uBound-m_included.getSum())) {
          node = m_included;
          return true;                        // Return signifying we found a feasible completion.
        }
      }
      continue;                               // Go to next node from stack
    }

    // If sumRemaining falls short of lower bound
    if (subsetSum + sumUnassigned < lowerBound) {
      s_tmpCount++;
      continue;                               // Go to next node from stack
    }

    // If all remaining elements fit in the set and form a feasible completion
    if (subsetSum + sumUnassigned <= m_uBound) {   // If all remaining elements fit in the set
      // In case the last element excluded is a duplicate of this new element, this case
      // would have happened already when the last element was included, so don't add again
      // TODO: Make sure this if check is a necessary step still
      if (m_excluded.empty() || (m_excluded.back() != currentElement)) {
        for (int i = newAssignedIdx; i < m_numElements; i++) {   // for each remaining element
          m_included.push_back(m_elements[i]);                  // include it in the set
        }
        if (!m_checkDominance ||
            !bp::isDominated(m_included.elements(), m_included.getSum(),
                             m_excluded.elements(),m_uBound-m_included.getSum())) {
          node = m_included;
          return true;     // return signifying we found a node.
        } else {
          continue;       // Node dominated, go to next node from stack
        }
      }
    }

    // ================================================================
    // If no termination conditions occur, push children onto the stack
    // ================================================================

    const uint64_t &newAssigned     = m_elements[newAssignedIdx];
    const uint64_t newSumUnassigned = sumUnassigned - newAssigned;

    // Push Exclusion first
//    uint64_t exclusionLowerBound = (subsetSum+newAssigned <= m_uBound) ?
//        std::max(lowerBound,subsetSum+newAssigned+1) : lowerBound;

    m_stack.push(IENode(newAssignedIdx, newSumUnassigned, subsetSum,
        numIncluded, numExcluded+1 , false, lowerBound));

    // Now inclusion but never allow inclusion to follow exclusion of duplicate value
    if (m_excluded.empty() || newAssigned != m_excluded.back()) {
      m_stack.push(IENode(newAssignedIdx, newSumUnassigned, subsetSum + newAssigned,
                          numIncluded+1, numExcluded , true, lowerBound));
    }

  } // end big while loop
  return false; // No more elements
}

size_t IECompletionGenerator::s_tmpCount= 0;
} /* namespace ss */
