/*
 * CreatePartitioningProblems.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: ethan
 */

#include "PackingUtils.hpp"
#include "ProgramOptionsUtils.hpp"
#include "CreateProblemUtils.hpp"
#include "ss/Schroeppel_Shamir.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <numeric> // std::accumulate

using namespace boost;
namespace po = boost::program_options;

using std::cout;
using std::endl;


double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

uint64_t computeSum(const std::vector<uint64_t> &S, const uint64_t mask) {

  uint64_t sum = 0;
  for (size_t i=0;i<S.size();i++) {  // For each element
    if ((mask >> i) & 1) {    // If bit is set in mask
       sum += S[i];           // Add corresponding element
    }
  }
  return sum;
}

void generateRandom(std::vector<uint64_t> &S, const size_t N,
										const uint64_t minValue, const uint64_t maxValue, const uint32_t seed) {
  // Setup random number generator
  typedef boost::mt19937 RNGType;
  RNGType randomGenerator;
  boost::uniform_int<int64_t> range( minValue, maxValue );
  boost::variate_generator< RNGType, boost::uniform_int<int64_t> >
                boostRandom(randomGenerator, range);
  boostRandom.engine().seed(seed);

  for (size_t i=0;i<N;i++) {
    S.push_back(boostRandom());       // Generate random elements
  }

  std::sort(S.begin(),S.end(),std::greater<int64_t>());   // sort elements descending

}
int main(int argc, char *argv[]) {
	std::cout.imbue(std::locale(""));
	if (argc != 2) {
		cout << "\n\n   Usage: " << argv[0] << " [N]" << endl << endl;
		exit(0);
	}
	const size_t  N           = atoi(argv[1]);   // The number of elements
  uint32_t seed = 0;
  const int     capBase     = 2;
  const int     capExponent = 20;

  const size_t NUM_SAMPLES = 5000000;

  const uint64_t minValue = 1;
  const uint64_t maxValue = ((int64_t) pow(capBase,capExponent))-1;  // maxValue = b^e - 1

//  cout << endl
//       << "Max Value: " << capBase << "^" << capExponent << " - 1 = " <<maxValue << endl
//       << "# Items  : " << N << endl << endl;


  std::vector<uint64_t> S;  // Vector to store elements
  generateRandom(S, N, minValue, maxValue, seed); 				// Generate random numbers

  uint64_t NUM_SUBSETS= ((uint64_t) 1) << N;														// 2^N subsets

  uint64_t totalSum = std::accumulate(S.begin(),S.end(),(uint64_t) 0);
  uint64_t mu = totalSum / N;



  // Compute sample variance
  // For k samples
  //
  // s^2 = 1/(k-1) x ([sum from i=1 to k]: (x_i - mu)^2)

  std::vector<uint64_t> masks;
  generateRandom(masks,NUM_SAMPLES,0,NUM_SUBSETS-1,seed);

//  uint64_t variance = 0;
//
//  for (size_t i=0;i<masks.size();i++) {
//  	uint64_t mask = masks[i];
//    uint64_t subsetSum  = computeSum(S,mask);
//    uint64_t sample = (subsetSum >= mu) ? subsetSum-mu : mu - subsetSum;
//    sample *= sample;
//    variance += sample;
//  }
//  variance /= (NUM_SAMPLES-1);
//  double std = sqrt(variance);
//  cout << "N: "        << std::setw(2)  << N
//       << "  Mu: "     << std::setw(10) << mu
//       << "  Var: "    << std::setw(12) << variance
//       << "  STD: "    << std::setw(10) << std
//  		 << "  Var/Mu: " << std::setw(2) << variance / mu << endl;
//
//  uint64_t avg = mu * N / 2;
//  uint64_t DIFF = avg / 10;
//  double z = (double) DIFF / (double) std;
//  cout << "Z: " << z << endl;
//  cout << "phi(z) = " << phi(z) -.5 << endl;
//  cout << "Num Samples: " << NUM_SUBSETS<< endl;
//  cout << "Prediction: " << NUM_SUBSETS * (phi(z)-.5) << endl;


  const size_t NUM_K = 15;
  const size_t MIN_K = 3;
  const size_t NUM_BUCKETS = 100;
  size_t histogram[NUM_BUCKETS][NUM_K];
  uint64_t bucketSize = totalSum / NUM_BUCKETS;
  memset(histogram,0,sizeof(size_t) * NUM_BUCKETS * NUM_K);

  for (size_t mask=0;mask<NUM_SUBSETS;mask++) {
  	int numBits = __builtin_popcount(mask);
  	if ( numBits >= MIN_K && numBits < NUM_K + MIN_K) {
			uint64_t subsetSum  = computeSum(S,mask);
			histogram[subsetSum / bucketSize][numBits-MIN_K]++;
  	}
  }

  for (size_t i=0;i<NUM_BUCKETS;i++) {
  	cout << std::setw(2) << i << ": ";

  	for (size_t j = 0; j<NUM_K;j++) {
  		cout << std::setw(9) << histogram[i][j] << "  ";
  	}
  	cout << endl;
  }
//  vector<ss::SSSetNode> sets;               // The sets generated by Schroeppel and Shamir
//
//  cout << "Min : " << avg << endl
//  		 << "Max : " << avg+DIFF << endl
//  		 << "Size: " << S.size() << endl;
//  ss::generateSetsSS(&S[0], S.size(), avg, avg+DIFF, sets);
//  cout << "Actual   : " << sets.size() << endl;
}


