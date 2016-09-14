/*
 * CreatePartitioningProblems.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: ethan
 */

#include "../pack/PackingUtils.hpp"
#include "ProgramOptionsUtils.hpp"
#include "CreateProblemUtils.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

using namespace boost;
namespace po = boost::program_options;

using std::cout;
using std::endl;

const string DIR_NAME = "./dat/partitioning/";


// ============================================================================
//
// R i c h R a n d o m
//
// ============================================================================
#define RICH_INIT_SEED  13070                    			/* initial random seed, in decimal */
#define ACONST 25214903917                            /* constant multiplier */
#define C 11                                          /* additive constant */
#define MASK 281474976710655                          /* 2^{48}-1 */
uint64_t g_richSeed;
uint64_t richRandom(uint64_t a[], int n, const uint64_t maxValue)

{int i;                                               /* index into array */
 long long sum;                                       /* sum of all numbers */

 sum = 0;                                    					/* initialize sum of all numbers */
 for (i = 0; i < n; i++)                         			/* for each element of array */
   {g_richSeed = (ACONST * g_richSeed + C) & MASK;    /* next seed in random sequence */
     a[i] = g_richSeed % maxValue;              			/* random value from zero to maxValue-1 */
     sum += a[i];}                              			/* compute sum of all numbers */
 return (sum);}


// ============================================================================
//
// M a i n
//
// ============================================================================
int main(int argc, char *argv[]) {

  uint32_t seed = 0;
  int     capBase     = UNSET_INT;
  int     capExponent = UNSET_INT;
  int     N           = UNSET_INT;   // The number of elements
  int     numTrials    = UNSET_INT;


  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"         , "produce help message")
  ("seed,s"         , po::value< uint32_t>(&seed)    , "The random number generator seed.")
  ("capBase,b"      , po::value<int>(&capBase)       , "(required) The capacity base.     (cap is b^e)")
  ("capExponent,e"  , po::value<int>(&capExponent)   , "(required) The capacity exponent. (cap is b^e)")
  ("numElements,n"  , po::value< int >(&N)           , "(required) The number of elements per instance.")
  ("numTrials,z"    , po::value< int >(&numTrials)    , "(required) The number of files to generate.");


  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("capBase"    , capBase     ,UNSET_INT);
    optionRequired("capExponent", capExponent ,UNSET_INT);
    optionRequired("numElements", N           ,UNSET_INT);
    optionRequired("numTrials"   , numTrials    ,UNSET_INT);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return 0;
    }
  } catch(std::exception& e) {
    cout << helpOut.str() << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return 1;
  }

  mkdir(DIR_NAME.c_str());    // In case dat dir does not exist

  const uint64_t minValue = 1;
  const uint64_t maxValue = ((int64_t) pow(capBase,capExponent))-1;  // maxValue = b^e - 1

  cout << endl
       << "Max Value: " << capBase << "^" << capExponent << " - 1 = " <<maxValue << endl
       << "# Items  : " << N << endl
       << "# Files  : " << numTrials << endl << endl;

  // Setup random number generator
  typedef boost::mt19937 RNGType;
  RNGType randomGenerator;
  boost::uniform_int<int64_t> range( minValue, maxValue );
  boost::variate_generator< RNGType, boost::uniform_int<int64_t> >
                boostRandom(randomGenerator, range);
  boostRandom.engine().seed(seed);

  // Create file
  std::ostringstream filename;
  filename << DIR_NAME << capBase << "_" << capExponent << "_" << N << "_" << numTrials << ".np";
  std::ofstream outFile(filename.str().c_str());
  cout << "Writing: " << filename.str() << endl;
  uint64_t *S = new uint64_t[N];  // Array to store elements
  g_richSeed = seed;
  for (int i=0;i<numTrials;i++) {
  	richRandom(S,N,maxValue+1);
//    for (int j=0;j<N;j++) {       // Generate random elements
//    	//S[j] = boostRandom();
//    }
    std::sort(S,S+N,std::greater<int64_t>());   // sort elements descending

    outFile << i       ;                        // Write problem index

    for (int j=0;j<N;j++) {                     // Write elements to file
      outFile << " " << S[j];
    }
    outFile << endl;                            // Go to next line for next problem
  }

  outFile.close();                              // Close file
  delete [] S;                                  // Free memory
}


