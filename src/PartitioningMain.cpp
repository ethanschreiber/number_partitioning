/*
 * MultiwayNumberPartitioning
 *
 *  Created on: Jul 18, 2012
 *      Author: ethan
 */

#include "main/MainUtils.hpp"
#include "Utils.hpp"
#include "partition/PartitionUtils.hpp"
#include "utils/ProgramOptionsUtils.hpp"
#include "ss/CGA.hpp"
#include "ss/CKK.hpp"
#include "ss/DPPartition.hpp"
#include "ss/SubsetSum.hpp"
#include "ss/SchroeppelShamir.hpp"
#include "ss/HorowitzSahni.hpp"
#include "partition/Moffitt.hpp"
#include "partition/SNP.hpp"
#include "partition/RNP.hpp"
#include "partition/PartitionIterativeWeakening.hpp"
#include "partition/PartitionIterativeWeakeningLowCardinality.hpp"
#include "partition/BinarySearch.hpp"


#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>   // For ceil
#include <stdint.h>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <sstream>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

using namespace boost;
namespace po = boost::program_options;

// Specification
bool readCommandLine(int argc,char *argv[],partition::PartitionOptions &partitionOptions);

// ============================================================================
// main reads the command line and calls execute.
// ============================================================================

int main(int argc, char *argv[])
{
  std::cout.imbue(std::locale("")); // For printing 1,000 instead of 1000
  // ==========================================================================
  // Read from command line
  // ==========================================================================

  partition::PartitionOptions partitionOptions;
  PackingOptions  packingOptions;

  bool shouldExit = readCommandLine(argc,argv,partitionOptions);
  if (shouldExit) { return 0; }

  SolutionMethodPartition solutionMethod = getSolutionMethodPartition(partitionOptions.solutionMethodInt);    // Set the solution method from the int
  std::vector<partition::PartitionProblem> problems = getProblems(partitionOptions);				// Read the problems from the problem files

  // ==========================================================================
  // Open the output file stream
  // ==========================================================================
  string  outputFilename = getOutputFilename(partitionOptions.inputFilename, getSuffix(solutionMethod),packingOptions, 2);

  int startIdx = (partitionOptions.problemIdx == UNSET_INT)	? // If no startIdx specified
  		countLines(outputFilename,"problem") :									// Get startIdx by counting lines in output file
  		partitionOptions.problemIdx;														// Otherwise start at the specified Index

  std::ofstream outFile;		// The output file, make sure to close it later
  std::streambuf *buf;			// The stream buffer
  openOutputFile(buf,outFile,outputFilename,startIdx,partitionOptions.problemIdx, problems.size());
  std::ostream out(buf);	// Use a streambuf so we can append to file

  // ==========================================================================
  // Loop through partition problems
  // ==========================================================================

  ProblemStats stats;
  size_t endIdx = (partitionOptions.problemIdx == UNSET_INT) ? problems.size() : startIdx+1;

  for (size_t i=startIdx;i<endIdx;i++) {
    partition::PartitionProblem &problem = problems[i];
//    cout << endl << "S : "; std::copy(problem.S,problem.S + problem.N, std::ostream_iterator<int>(std::cout," "));
    uint64_t result;
    SimpleTimer timer;
    std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order

    if (solutionMethod == CGA) {
    	result = ss::executeCGA2(problem.S,problem.N,problem.sum);

    } else if (solutionMethod == CKK) {
      result = ss::executeCKK2(problem.S,problem.N,problem.sum);

    } else if (solutionMethod == HS) {
      result = ss::executeHS2(problem,stats);
    } else if (solutionMethod == SS) {

      result = ss::executeSS2(problem,stats);
    } else if (solutionMethod == DP) {
    	result = ss::executeDPPartition(problem,stats);
    } else {
      cout << "This should never happen, unknown solution method!" << endl;
      exit(0);
    }

    uint64_t perfect = (problem.sum+1) / 2;
    bool isPerfect = (perfect == result);
    if (solutionMethod == HS ||
    		solutionMethod == SS ||
    		solutionMethod == DP) {
    	out << problem.problemName <<  " " << timer.timeElapsed()  <<  " " <<  result << " " << isPerfect << " " << std::fixed << std::setprecision(2) << stats.residentMemory << endl;
    } else {
    	out << problem.problemName <<  " " << timer.timeElapsed()  <<  " " <<  result << " " << isPerfect << endl;
    }


  }

  outFile.close();
  return 0;
}

// ================================================================
// Read from command line. Returns true if the program should
// exit due to help being asked for or some error, false otherwise.
// ================================================================
bool readCommandLine(int argc, char *argv[], partition::PartitionOptions &partitionOptions) {

	std::ostringstream solutionMethodString;
	solutionMethodString << "(required) The solution method (" << partitionMethodsToString() << ")";

  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"           , "produce help message")
  ("method,m"         , po::value< int >    (&partitionOptions.solutionMethodInt) , solutionMethodString.str().c_str())
  ("file,f"           , po::value< string > (&partitionOptions.inputFilename)     , "(required) The input filename.")
  ("problem-index,x"  , po::value< int >    (&partitionOptions.problemIdx)        , "The problem index from the input filename. If this is left out, all problems in the file are solved");


  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("method"					, partitionOptions.solutionMethodInt, UNSET_INT);
    optionRequired("file"						, partitionOptions.inputFilename		, UNSET_STRING);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return true;
    }
  } catch(std::exception& e) {
    cout << helpOut.str()
         << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return true;
  }
  return false;
}


