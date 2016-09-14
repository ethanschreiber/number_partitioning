/*
 * PartitionUtils.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: ethan
 */

#include "PartitionUtils.hpp"


// *******************
// PartitioningProblem
// *******************
namespace partition {
  PartitionProblem::PartitionProblem(int N0, int K0, uint64_t S0[], uint64_t sum0, string problemName0 )
    : PackingProblem(N0,S0,sum0,problemName0), K(K0) {}


  PartitionProblem::PartitionProblem(const string &filename) {
    load(filename);
  }

  PartitionProblem::PartitionProblem(const PartitionProblem &problem)
  : PackingProblem(problem), K(problem.K) {}

  PartitionProblem::PartitionProblem(const PackingProblem &problem, int K0)
  : PackingProblem(problem), K(K0) {}

  void PartitionProblem::load(const string &filename) {
    std::ifstream inFile(filename);
    std::getline(inFile,problemName);   // Read problem Name

    inFile >> K >> N >> sum;            // [k] [N] [sum]
    readElements(inFile);

  }
} //end namespace
