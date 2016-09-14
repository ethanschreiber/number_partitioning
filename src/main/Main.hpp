/*
 * Main.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef MAIN_HPP_
#define MAIN_HPP_

#include "CLI.hpp"
#include "MainUtils.hpp"
#include "../Utils.hpp"

/**
 * Main class to run
 */
class Main {
protected:
	virtual SolutionMethodKPartition getSolutionMethod() = 0;
	virtual void init(partition::PartitionOptions &partitionOptions) { }
	virtual uint64_t executePartition(const partition::PartitionProblem &problem, const partition::PartitionOptions &partitionOptions,
																		const PackingOptions &packingOptions, ProblemStats &stats) = 0;

	// Use this to both update values if necesary and return a string which is to be
	// printed to the output file. THis is overridden for CIW and CIW_LC
	virtual string processResults(partition::PartitionOptions &partitionOptions, const ProblemStats stats,
			partition::PartitionProblem &problem, const double timeElapsed, const uint64_t result, const int idx) {
		std::ostringstream out;
		out << problem.problemName <<  " " << timeElapsed <<  " " <<  result << endl;
		return out.str();
	}
public :
	virtual ~Main() {}




	int execute(int argc, char *argv[]) {
	  // ==========================================================================
	  // Read from command line
	  // ==========================================================================

	  partition::PartitionOptions partitionOptions;
	  PackingOptions packingOptions;
	  bool shouldExit = cli::readCommandLine(argc,argv,partitionOptions, packingOptions, getSolutionMethod());
	  if (shouldExit) { return 0; }

	  std::vector<partition::PartitionProblem> problems = getProblems(partitionOptions);				// Read the problems from the problem files
	  init(partitionOptions);
	  // ==========================================================================
	  // Open the output file stream
	  // ==========================================================================

	  string  outputFilename = getOutputFilename(partitionOptions.inputFilename, getSuffix(getSolutionMethod()),packingOptions, partitionOptions.inputK);

	  int startIdx = (partitionOptions.problemIdx == UNSET_INT)	? // If no startIdx specified
	  		countLines(outputFilename,"problem") :									// Get startIdx by counting lines in output file
	  		partitionOptions.problemIdx;														// Otherwise start at the specified Index
	  size_t endIdx = (partitionOptions.problemIdx == UNSET_INT) ? problems.size() : startIdx+1;

	  std::ofstream outFile;		// The output file, make sure to close it later
	  std::streambuf *buf;			// The stream buffer
	  openOutputFile(buf,outFile,outputFilename,startIdx,partitionOptions.problemIdx, problems.size());
	  std::ostream out(buf);	// Use a streambuf so we can append to file

	  // ==========================================================================
	  // Loop through partition problems
	  // ==========================================================================
	  ProblemStats stats;

	  for (size_t i=startIdx;i<endIdx;i++) {
	    partition::PartitionProblem &problem = problems[i];
	    SimpleTimer timer;

	    uint64_t result = executePartition(problem, partitionOptions, packingOptions, stats);

	    out << processResults(partitionOptions, stats, problem, timer.timeElapsed(), result, i);
	    out.flush();
	  }

	  outFile.close();


	  return 0;

	}
};

#endif /* MAIN_HPP_ */
