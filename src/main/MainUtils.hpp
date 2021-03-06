/*
 * MainUtils.hpp
 *
 *  Created on: Mar 14, 2014
 *      Author: ethan
 */

#ifndef MAINUTILS_HPP_
#define MAINUTILS_HPP_

#include "../partition/PartitionUtils.hpp"
#include <string>
#include <vector>

using std::string;
using std::vector;

// =============================================
// Function Specifications for utility functions
// =============================================
std::vector<partition::PartitionProblem> getProblems(const partition::PartitionOptions &partitionOptions);
std::vector<uint64_t> getSolutionValues(const string filename);
void openOutputFile(std::streambuf *& buf, std::ofstream &outFile, const string outputFilename,const int startIdx,
									  const int problemIdx, const size_t numProblems);

#endif /* MAINUTILS_HPP_ */
