/*
 * MultiwayNumberPartitioning
 *
 *  Created on: Jul 18, 2012
 *      Author: ethan
 */

#include "SubsetSumMain.hpp"

#include <deque>
#include <boost/program_options.hpp>
#include "utils/ProgramOptionsUtils.hpp"
using namespace boost;
namespace po = boost::program_options;


// ============================================================================
// main reads the command line and calls execute.
// ============================================================================
int main(int argc, char *argv[])
{
  string inputFilename = UNSET_STRING;

  double ubPercent = UNSET_DOUBLE;
  double lbPercent = UNSET_DOUBLE;

  bool isDynamicProgramming   = false;
  bool isHorowitzSahni        = false;
  bool isInclusionExclusion   = false;
  bool isOrderedPowerSet      = false;
  bool isSchroeppelShamir     = false;

  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"               , "produce help message")
  ("file,f"               , po::value< string >(&inputFilename), "(required) The input filename.")
  ("lower-bound-percent,l", po::value< double >(&lbPercent), "The lower bound for the binary search.")
  ("upper-bound-percent,u", po::value< double >(&ubPercent), "The upper bound for the binary search.")
  ("dp,d"                 , po::value(&isDynamicProgramming)->zero_tokens(),"Should we use dynamic programming?")
  ("hs,h"                 , po::value(&isHorowitzSahni     )->zero_tokens(),"Should we use Horowitz and Sahni?")
  ("ie,i"                 , po::value(&isInclusionExclusion)->zero_tokens(),"Should we use Inclusion/Exclusion?")
  ("ops,o"                , po::value(&isOrderedPowerSet   )->zero_tokens(),"Should we use ordered power set?")
  ("ss,s"                 , po::value(&isSchroeppelShamir  )->zero_tokens(),"Should we use Schroeppel and Shamir?");
  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("file",inputFilename,UNSET_STRING);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return 0;
    }
  } catch(std::exception& e) {
    cout << helpOut.str()
         << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return 1;
  }

  if (ubPercent == UNSET_DOUBLE || ubPercent > 1 || ubPercent < 0) {
    ubPercent = 1.0;
  }

  if (lbPercent == UNSET_DOUBLE || lbPercent > 1 || lbPercent < 0) {
    lbPercent = 0;
  }


  vector<string> inputFilenames = getInputFilenames(inputFilename);
  string outputFilename = getOutputFilename(inputFilename,lbPercent,ubPercent);
  ofstream out(outputFilename);
  cout  << endl
        << "Input Filename : " << inputFilename << endl
        << "Output Filename: " << outputFilename << endl << endl;



  cout << "Num Files: " << inputFilenames.size() << endl;
  for (vector<string>::iterator it = inputFilenames.begin();it != inputFilenames.end(); it++) {
    PackingProblem P(*it);
    std::sort( P.S, P.S+P.N ,std::less<uint64_t>());  // Sort in ascending order

    out << P.problemName <<  " ";


    vector<size_t> numSubsets;
    uint64_t upperBound = (uint64_t) std::ceil((double) P.sum * ubPercent);
    uint64_t lowerBound = (uint64_t) std::floor((double) P.sum * lbPercent);

    cout << P.problemName << ": (Sum = " << P.sum
         << ") (LB = " << lowerBound
         << ") (UB = " << upperBound << ")...";
    cout.flush();
    // -------------------
    // Inclusion Exclusion
    // -------------------
    if (isInclusionExclusion) {
      cout << "Inclusion/Exclusion  ";
      cout.flush();
      SimpleTimer timer;
      vector<ss::SetNodeBitset> sets;
      size_t count = ss::incExc(P.S,P.N,P.sum,lowerBound,upperBound,sets);
      numSubsets.push_back(count);
      cout << "IE: " << count << endl;
      out << "(IE " << timer.timeElapsed() << ") " ;
    } else {
      out << "-1 ";
    }

    // ---------------------
    // Schroeppel and Shamir
    // ---------------------
    if (isSchroeppelShamir) {
      cout << "Schroeppel and Shamir  ";
      cout.flush();
      SimpleTimer timer;
      vector<ss::SetNodeBitset> sets;
      size_t count = ss::ESSSets(P.S,P.N,lowerBound,upperBound,sets);
      numSubsets.push_back(count);
      out << "(SS " << timer.timeElapsed() << ") ";
    } else {
      out << "-1 ";
    }

    // ------------------
    // Horowitz and Sahni
    // ------------------
    if (isHorowitzSahni) {
      cout << "Horowitz and Sahni  ";
      cout.flush();
      SimpleTimer timer;
      vector<ss::SetNodeBitset> sets;
      size_t count = ss::EHS(P.S,P.N,lowerBound,upperBound,sets);
      numSubsets.push_back(count);
      out << "(HS " << timer.timeElapsed() << ") ";
    } else {
      out << "-1 ";
    }
    // -----------------
    // Ordered Power Set
    // -----------------
    if (isOrderedPowerSet) {
      cout << "Ordered Power Set  ";
      cout.flush();
      SimpleTimer timer;
      std::deque<ss::SetNodeBitset> sets;
      size_t count = ss::generateSetsOPS(P.S,P.N,lowerBound,upperBound,sets);
      numSubsets.push_back(count);
      cout << "OPS: " << count << endl;
      out << "(OPS " << timer.timeElapsed() << ") ";
    } else {
      out << "-1 ";
    }

    // -------------------
    // Dynamic Programming
    // -------------------
//    if (isDynamicProgramming) {
//      cout << "Dynamic Programming  ";
//      cout.flush();
//      SimpleTimer timer;
//      vector<ss::Bitset> sets;
//      ss::DPMatrix matrix(P.S,P.N,P.S[0],upperBound);
//      size_t count = matrix.generateSubsets(sets,lowerBound,upperBound);
//
//      numSubsets.push_back(count);
//      out << "(DP " << timer.timeElapsed() << ") ";
//    } else {
//      out << "-1 ";
//    }
//    cout << "DONE" << endl;
//

    // -----------------
    // Now Check Results
    // -----------------
    size_t numSets = numSubsets[0];
    for (size_t i=0;i<numSubsets.size();i++) {
      if (numSubsets[i] != numSets) {
        cout << "ERROR, inconsistent results: " << numSets << " != " << numSubsets[i] << endl;
        //exit(0);
      }
    }
    out << numSets << endl;

  }
  out.close();
  return 0;
}
