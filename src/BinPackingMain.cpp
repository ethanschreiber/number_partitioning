/*
 * BinPackingRichMainSingle.cpp
 *
 *  Created on: Jul 18, 2012
 *      Author: ethan
 */

#include "pack/PackingUtils.hpp"
#include "utils/ProgramOptionsUtils.hpp"
#include "pack/BinCompletion.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

// Specification
bool readCommandLine(int argc,char *argv[],PackingOptions &packingOptions);


// ============================================================================
// main reads the command line and calls execute.
// ============================================================================
int main(int argc, char *argv[])
{
  // ===============
  // Process Options
  // ===============

  PackingOptions packingOptions;

  bool shouldExit = readCommandLine(argc,argv,packingOptions);
  if (shouldExit) { return 0; }

  SolutionMethodKPartition solutionMethod = getSolutionMethodKPartition(packingOptions.solutionMethodInt);
  string outputFilename = getOutputFilename(packingOptions.inputFilename,getSuffix(solutionMethod), packingOptions,UNSET_INT);
  vector<string> inputFilenames = getInputFilenames(packingOptions.inputFilename);

  // ===========================================================
  // Check to see if some of the file already has been completed
  // ===========================================================
  int startIdx =
   		(outputFilename.find("_all_") == string::npos) ? // If for single file
 			0 :																							 // compute it even if its been done
   		countLines(outputFilename,"problem");						// Otherwise pickup if we are in the middle of a dir

  // Open output file in append mode
  std::ofstream outFile(outputFilename.c_str(), std::ios_base::app | std::ios_base::out);

  {
    std::ostringstream os;
    os << "Output Filename: " << outputFilename << " ";
    if (startIdx >= inputFilenames.size()) {
      os << "[File already complete!]";
    } else if (startIdx > 0) {
      os << "[Appending from problem " << startIdx << ".]";
    }
    cout << os.str() << endl;
  }
  // ========================
  // Loop through input files
  // ========================
  ProblemStats stats;
  for (size_t i=startIdx;i<inputFilenames.size();i++) {
    string filename = inputFilenames[i];

    BinPackingProblem *problem = new BinPackingProblem(filename);


    string problemName(problem->problemName);
    if (packingOptions.isVerbose) {
      cout << "------- " << problemName << " -------" << endl;
    }
    if (solutionMethod == BSBC) {
      executeBP(*problem,stats,packingOptions);
    } else if (solutionMethod == BSBCP) {
      executeBelovBCP(*problem,stats);
    } else {
      executeTest(*problem,stats,packingOptions);
    }

    outFile << problemName  << ","
            << stats.time           << ","
            << stats.numNodes       << ","
            << stats.numBins        << ","
            << stats.sum            << ","
            << stats.lowerbound     << endl;

    if (packingOptions.isVerbose) {
      cout
         << "Time           : " << stats.time           << endl
         << "Node Count     : " << stats.numNodes       << endl
         << "Bin Count      : " << stats.numBins        << endl
         << "Sum            : " << stats.sum            << endl
         << "Lower Bound    : " << stats.lowerbound     << endl
         << "# Perfect Pairs: " << stats.numPerfectPairs << endl;
    }

      delete problem;
  }

  outFile.close();

}


// ================================================================
// Read from command line. Returns true if the program should
// exit due to help being asked for or some error, false otherwise.
// ================================================================
bool readCommandLine(int argc, char *argv[],PackingOptions &packingOptions) {


  po::options_description commandOptions("Options");
  commandOptions.add_options()
		("help,h"                                                                         , "produce help message"                                                        )
		("method,m"       , po::value <int   > (&packingOptions.solutionMethodInt)		    , "(required) The solution method (0=Bin-Completion, 2=Branch-and-Cut-and-Price")
		("buffer-size,b"  , po::value <int   > (&packingOptions.bufferSize)               , "(Default=50) The size of the buffer for inclusion/exclusion search"          )
		("file,f"         , po::value <string> (&packingOptions.inputFilename)						, "(required) The input filename."                                              )
		("lds,d"          , po::value(&packingOptions.useLDS)->zero_tokens()							, "Should we use limited discrepency search? (flag)"                            )
		("printsol,p"     , po::value(&packingOptions.printSolution)->zero_tokens()			  , "Should we print the solution? (flag)"                                        )
		("verbose,v"      , po::value(&packingOptions.isVerbose)->zero_tokens()					  , "Should we display verbose output? (flag)"                                    )
		("remove-pairs,r" , po::value(&packingOptions.removePairs)->zero_tokens()				  , "Should we remove pairs during preprocess stage? (flag)"                      )
		("ss,s"           , po::value(&packingOptions.useSchroeppelShamir)->zero_tokens() , "Should we use Schroeppel and Shamir to generate completions? (flag)"         )
		("classic-sort,c" , po::value(&packingOptions.classicSort)->zero_tokens()				  , "Should we use the classic korf bin completion sort?"                         );

  std::ostringstream helpOut;
  helpOut << "\nUsage: " << argv[0] << " [options]\n\n" << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("method",packingOptions.solutionMethodInt,UNSET_INT);
    optionRequired("file"	 ,packingOptions.inputFilename,UNSET_STRING);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return true;
    }
  } catch(std::exception& e) {
    cout << helpOut.str() << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return true;
  }
  return false;
}
