/*
 * Moffitt.hpp
 *
 *  Created on: Aug 20, 2013
 *      Author: ethan
 */

#ifndef CGA_MULTIWAY_HPP_
#define CGA_MULTIWAY_HPP_

#define PAPER_DEBUG 0

#include "PartitionUtils.hpp"

// CGA Multiway.

namespace partition {

class CGA_MW {
private :
	uint64_t m_best;
	int m_n;
	int m_k;
	uint64_t *m_S;
	uint64_t m_sum;
protected :
	uint64_t search(int depth, uint64_t *&subsets, uint64_t sum);

	void swapSort(uint64_t *&subsets,int idx);
public :

	CGA_MW(const PartitionProblem &problem);
	uint64_t search();
	~CGA_MW();

};

uint64_t executeCGA_MW(const PartitionProblem &problem, ProblemStats &stats);

} // End Namespace

#endif /* MOFFITT_HPP_ */
