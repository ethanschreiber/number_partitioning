/*
 * BinarySearch.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: ethan
 */

#include "BinarySearch.hpp"
#include "../ss/extended/Extended_Schroeppel_Shamir.hpp"
namespace binary_search {


uint64_t BinarySearch::search(BinPackingProblem &problem, uint64_t minCapacity,uint64_t maxCapacity,
                     	 	 	 	 	 	 	 int numPartitions,  uint64_t maxElement,
                     const PackingOptions &packingOptions) {

  std::cout.imbue(std::locale("")); // For printing 1,000 instead of 1000
  ProblemStats stats;

//  cout << "Min: " << minCapacity << endl
//       << "Max: " << maxCapacity << endl << endl;
  // ----------------------------------------------------------------------------------
  // If the maxElement is greater than 2^32, use schroeppel and shamir to only search
  // possible subset sums. The idea is the possible values is sparse for large numbers
  // ----------------------------------------------------------------------------------
  //cout << endl;

  if (packingOptions.useSchroeppelShamirBinarySearch) {

//    cout << "Using SS" << endl;
    vector<uint64_t> SSSets;
    int minIdx = 0;

    std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order
    while (maxCapacity - minCapacity > 10000) {
      cout << minCapacity << " " << maxCapacity << " " << std::setw(15) << maxCapacity - minCapacity << endl;

      uint64_t midCapacity = (maxCapacity + minCapacity) / 2;
      problem.capacity    = midCapacity;     // Set capacity


      execute(problem,stats,packingOptions);

      if (stats.numBins > numPartitions) {
        minCapacity = midCapacity + 1;
      } else {
        maxCapacity = midCapacity;
        //maxCapacity = stats.maxCapUsed;
      }
    }

//    cout << "Doing SS" << endl;
    // Find the min and max capacity
    std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in ascending order
    int maxIdx = ss::ESSSums(problem.S,problem.N,minCapacity,maxCapacity,SSSets)-1;
    std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order

    // ETHAN: This is a possible performance bottleneck!
    std::sort( SSSets.begin(), SSSets.end(),std::less<uint64_t>());        // Sort results in ascending order



    while (maxIdx > minIdx) {
      //cout << SSSets[minIdx] << " " << SSSets[maxIdx] << " " << SSSets[maxIdx] - SSSets[minIdx] << endl;
      int midIdx       = (maxIdx + minIdx) / 2;
      problem.capacity = SSSets[midIdx];

      execute(problem,stats,packingOptions);

      if (stats.numBins > numPartitions) {
        minIdx = midIdx + 1;
      } else {
        maxIdx = midIdx;
        //maxCapacity = stats.maxCapUsed;
      }
    }

    return SSSets[minIdx];


    // ----------------------------------------------------------------------------------
    // For small numbers with maxElement < 2^32, don't use Schroeppel and Shamir, just
    // do a direct binary search
    // ----------------------------------------------------------------------------------
  } else {

    std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order
    while (maxCapacity > minCapacity) {
//      cout << minCapacity << " " << maxCapacity << " " << std::setw(15) << maxCapacity - minCapacity << endl;

      uint64_t midCapacity = (maxCapacity + minCapacity) / 2;
      problem.capacity    = midCapacity;     // Set capacity

      execute(problem,stats,packingOptions);

      if (stats.numBins > numPartitions) {
        minCapacity = midCapacity + 1;
      } else {
        maxCapacity = midCapacity;
      }

//      // TEMP
//      uint64_t CMin = problem.sum - (numPartitions-1) * (maxCapacity);
//      vector<uint64_t> SSSets;
//      std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in ascending order
//      int maxIdx = ss::generateSumsSS(problem.S,problem.N,CMin,maxCapacity,SSSets)-1;
//      std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order
//
//      cout << "# SS Sets: " << SSSets.size() << endl
//           << "      Min: " << CMin << endl
//           << "      Max: " << maxCapacity << endl << endl;
//      // END TEMP

    }
    return minCapacity;
  }
}

} // end namespace
