/*
 * CreateProblems.cpp
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


const string DAT_BASE = "./dat/packing/";

int main(int argc, char *argv[]) {

  uint32_t seed = 0;
  int     numFiles    = UNSET_INT;
  int     N           = UNSET_INT;   // The number of elements
  int     capBase     = UNSET_INT;
  int     capExponent = UNSET_INT;
  double  maxFraction = UNSET_DOUBLE;

  uint64_t capacity;  // capacity or partition
  const uint64_t minValue = 1;
  uint64_t maxValue;

  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"       , "produce help message")
  ("seed,s"       , po::value< uint32_t>(&seed)    , "The random number generator seed.")
  ("numFiles,z"   , po::value< int >(&numFiles)    , "(required) The number of files to generate.")
  ("numElements,n", po::value< int >(&N)           , "(required) The number of elements per instance.")
  ("capBase,b"    , po::value<int>(&capBase)       , "(required) The capacity base.     (cap is b^e)")
  ("capExponent,e", po::value<int>(&capExponent)   , "(required) The capacity exponent. (cap is b^e)")
  ("maxFraction,f", po::value<double>(&maxFraction), "(required) The maximum value fraction (max is f * cap).");

  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("numFiles"   , numFiles    ,UNSET_INT);
    optionRequired("numElements", N           ,UNSET_INT);
    optionRequired("capBase"    , capBase     ,UNSET_INT);
    optionRequired("capExponent", capExponent ,UNSET_INT);
    optionRequired("maxFraction", maxFraction ,UNSET_DOUBLE);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return 0;
    }
  } catch(std::exception& e) {
    cout << helpOut.str() << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return 1;
  }


  capacity = (uint64_t) pow(capBase,capExponent);  // capacity = b^e
  maxValue = maxFraction * capacity;              // maxValue = b^e * f

  std::ostringstream dirName;

  cout << endl
       << "# Files  : " << numFiles << endl
       << "# Items  : " << N << endl
       << "Capacity : " << capacity << endl
       << "Max Value: " << maxFraction << " * " << capBase << "^" << capExponent << " = " <<maxValue << endl << endl;

  dirName << DAT_BASE << N << "/"
          << capBase << "_" << capExponent << "/"
          << maxFraction << "/";

  cout << dirName.str() << " - ";
  mkdir(dirName.str().c_str());


  // Setup random number generator
  typedef boost::mt19937 RNGType;
  RNGType randomGenerator;
  boost::uniform_int<uint64_t> range( minValue, maxValue  );
  boost::variate_generator< RNGType, boost::uniform_int<uint64_t> >
                boostRandom(randomGenerator, range);
  boostRandom.engine().seed(seed);

  size_t totalCount = 0;
  for (int i=0;i<numFiles;i++) {
     std::ostringstream problemName;
     std::ostringstream filename;
     problemName << "problem" << std::setfill('0') << std::setw(3) << i;
     filename << dirName.str() << problemName.str() << ".bp";

     uint64_t S[MAXN];

     uint64_t sum;
     int lowerBound;
     int upperBound;
     //uint64_t waste;
     //double wasteRatio;

     // -------------------------
     // Fill in the random values
     // -------------------------
     do {
       sum = 0;
       for (int j=0;j<N;j++) {
         S[j] = boostRandom();
         sum += S[j];
       }


         std::sort(S,S+N,std::greater<uint64_t>());

       int bfdDummy[MAXN];
       lowerBound = L2(S,N,capacity,sum);
       upperBound = BFD(S,bfdDummy,N,capacity);
       //waste = lowerBound * capacity - sum;
       //wasteRatio = (double) waste / (double) capacity;

       totalCount++;
     } while (upperBound == lowerBound); //|| wasteRatio > .00001); // TODO: Make this a parameter


     if (i % 5 == 0) {
       cout << i << " " ;
     }
     cout.flush();
     std::ofstream outFile(filename.str().c_str());

     outFile << problemName.str() << endl
             << capacity << " " << N << " " << sum << endl;



     for (int j=0;j<N;j++) {
       outFile << S[j] << endl;
     }

     outFile << endl;
     outFile.close();
   }
   cout << " DONE" << endl;
   cout << "Percent Trivial: " << std::fixed << std::setprecision(2) << ((double) (totalCount-numFiles) / (double) totalCount) * 100 << endl;

}
