/*
 * MultiwayNumberPartitioning
 *
 *  Created on: Jul 18, 2012
 *      Author: ethan
 */

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include "globals.hpp"
#include "main/MainUtils.hpp"
#include "pack/PackingUtils.hpp"
#include "partition/BinarySearch.hpp"
#include "partition/Moffitt.hpp"
#include "partition/PartitionIterativeWeakening.hpp"
#include "partition/PartitionIterativeWeakeningLowCardinality.hpp"
#include "partition/PartitionUtils.hpp"
#include "partition/RNP.hpp"
#include "partition/SNP.hpp"
#include "ss/CGA.hpp"
#include "utils/ProgramOptionsUtils.hpp"
#include "Utils.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

using namespace boost;
namespace po = boost::program_options;

// Specification
bool readCommandLine(int argc,char *argv[],partition::PartitionOptions &partitionOptions,PackingOptions &packingOptions);

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

  bool shouldExit = readCommandLine(argc,argv,partitionOptions,packingOptions);
  if (shouldExit) { return 0; }

  SolutionMethodKPartition solutionMethod = getSolutionMethodKPartition(partitionOptions.solutionMethodInt);    // Set the solution method from the int
  std::vector<partition::PartitionProblem> problems = getProblems(partitionOptions);				// Read the problems from the problem files

  // ==========================================================================
  // Open the output file stream
  // ==========================================================================
  string  outputFilename = getOutputFilename(partitionOptions.inputFilename, getSuffix(solutionMethod),packingOptions, partitionOptions.inputK);

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

    uint64_t result;
    SimpleTimer timer;

    //---------------------------------------------------------------------------------------
    if ((solutionMethod == BSBC) || (solutionMethod == BSBCP)) {
    //---------------------------------------------------------------------------------------

//      std::sort( problem.S, problem.S+problem.N ,std::greater<uint64_t>());  // Sort input in descending order
//
//      // Find the max input element
//      uint64_t maxElement = *std::max_element(problem.S,problem.S + problem.N);
//
//      // If the minCapacity was not passed from the command line
//      if (partitionOptions.minCapacity == UNSET_UINT64_T) {
//        // Can't be smaller than ([the sum of the elements] / [# partitions]) and also can't be smaller than maxElement
//        partitionOptions.minCapacity = std::max((uint64_t) ceil((double) problem.sum / (double) problem.K),maxElement);
//
//        // Also can't be smaller than the kth + k+1st element
//        partitionOptions.minCapacity = std::max(partitionOptions.minCapacity,problem.S[problem.K] + problem.S[problem.K]);
//      }
//
//      if (partitionOptions.maxCapacity == UNSET_UINT64_T) {
//        partitionOptions.maxCapacity = kk(problem.S,problem.N,problem.K,problem.sum);
//      }
//
//      BinPackingProblem binPackingProblem(problem,partitionOptions.maxCapacity);
//
//      result = binary_search::binarySearch(binPackingProblem,partitionOptions.minCapacity,partitionOptions.maxCapacity+1,
//      																		 problem.K,solutionMethod,maxElement,packingOptions);

    // --------------------------------------------------
    } else if (solutionMethod == RNP) {
    // --------------------------------------------------
      if (problem.N > MAXN) {
        cout << endl
             << "ERROR Problem Size: " << problem.N << endl
             << "  MAX Problem Size: " << MAXN << endl << endl;
        exit(0);
      }


    	result = executeRNP(problem,stats);
    } else if (solutionMethod == MOFFITT) {
      result = partition::executeMoffitt(problem,stats);
    } else if (solutionMethod == CIW) {
    	result = partition::executeCIW(problem, stats,partitionOptions.numSets);

      // DEBUG START

    	// update numSetscout << "KPART" << endl;
      if (i == 0) {
        partitionOptions.numSets = stats.firstCount;
//        cout << "New Num Sets: " << partitionOptions.numSets << endl;
      } else{
        partitionOptions.numSets = std::max(partitionOptions.numSets,stats.firstCount);
//        cout << "New Num Sets: " << partitionOptions.numSets << endl;
      }


      // DEBUG END
    } else if (solutionMethod == SNP) {
      if (problem.N > MAXN) {
        cout << endl
             << "ERROR Problem Size: " << problem.N << endl
             << "  MAX Problem Size: " << MAXN << endl << endl;
        exit(0);
      }

      result = snp::executeSNP(problem, stats);
    } else if (solutionMethod == CIW_LC) {
    	result = partition::executeCIWLowCardinality(problem,stats,partitionOptions.numSets);
    } else {
      cout << "This should never happen, unknown solution method!" << endl;
      exit(0);
    }

    if (solutionMethod == CIW || solutionMethod == CIW_LC) {
    	out << problem.problemName <<  " "
    			<< timer.timeElapsed()  <<  " "
    			<< result << " "
    			<< stats.firstCount << " "
    			<< stats.ssCalls << " "
    			<< stats.ssTime << " "
    			<< std::fixed << std::setprecision(2) << stats.residentMemory << endl;

/*
    	cout << problem.problemName <<  " "
    			<< timer.timeElapsed()  <<  " "
    			<< result << " "
    			<< stats.firstCount << " "
    			<< stats.ssCalls << " "
    			<< stats.ssTime << " "
    			<< std::fixed << std::setprecision(2) << stats.residentMemory << endl << "---------------------------" << endl;
*/
    } else {
    	out << problem.problemName <<  " " << timer.timeElapsed()  <<  " " <<  result << endl;
    }

//    cout << "Answer: " << result << endl;
//    cout << "Time  : " << timer.timeElapsed() << endl;

#if OLD_TECHNOLOGY
    minCapacity = UNSET_UINT64_T;    // For next iteration, make sure minCapacity is UNSET
    maxCapacity = UNSET_UINT64_T;    // For next iteration, make sure maxCapacity is UNSET
#endif
  }

  outFile.close();


//  cout << "Count: " << ss::IECompletionGenerator::s_tmpCount << endl;
  return 0;
}

// ================================================================
// Read from command line. Returns true if the program should
// exit due to help being asked for or some error, false otherwise.
// ================================================================
bool readCommandLine(int argc, char *argv[], partition::PartitionOptions &partitionOptions,PackingOptions &packingOptions) {

	std::ostringstream solutionMethodString;
	solutionMethodString << "(required) The solution method (" << kPartitionMethodsToString() << ")";

  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"           , "produce help message")
  ("num-partitions,k" , po::value< int >    (&partitionOptions.inputK)            , "The number of partitions. (required)")
  ("method,m"         , po::value< int >    (&partitionOptions.solutionMethodInt) , solutionMethodString.str().c_str())
  ("file,f"           , po::value< string > (&partitionOptions.inputFilename)     , "(required) The input filename.")
  ("problem-index,x"  , po::value< int >    (&partitionOptions.problemIdx)        , "The problem index from the input filename. If this is left out, all problems in the file are solved")
  ("num-cached-sets,n", po::value< size_t > (&partitionOptions.numSets)           ,"(DEFAULT=20000) The number of cached sets to generate for each schroeppel and shamir run. THis is used for -m 7 and -m 8.");

  // For Binary Search BC or BCP
  commandOptions.add_options()
  ("lower-bound,l"    , po::value< uint64_t >(&partitionOptions.minCapacity)   , "The lower bound for the binary search.")
  ("upper-bound,u"    , po::value< uint64_t >(&partitionOptions.maxCapacity)   , "The upper bound for the binary search.")
  ("buffer-size,b"    , po::value< int >(&packingOptions.bufferSize)            , "(Default=50) The size of the buffer for inclusion/exclusion search")
  ("lds,d"            , po::value(&packingOptions.useLDS)                          ->zero_tokens(),"Should we use limited discrepency search (flag, include or don't include.")
  ("ss,s"             , po::value(&packingOptions.useSchroeppelShamir)             ->zero_tokens(),"Should we use Schroeppel and Shamir to generate completions? (flag)")
  ("ssbs,t"           , po::value(&packingOptions.useSchroeppelShamirBinarySearch) ->zero_tokens(),"Should we use Schroeppel and Shamir to generate possible values for the binary search? (flag)")
  ("printsol,p"       , po::value(&packingOptions.printSolution)                   ->zero_tokens(),"Should we print the solution? (flag)")
  ("verbose,v"        , po::value(&packingOptions.isVerbose)                       ->zero_tokens(),"Should we display verbose output? (flag)")
  ("classic-sort,c"   , po::value(&packingOptions.classicSort)->zero_tokens(),"Should we use the classic korf bin completion sort?") ;

  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("num-partitions"	, partitionOptions.inputK						, UNSET_INT);
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


