/*
 * cli.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#include "CLI.hpp"


typedef po::options_description Options;



namespace cli {


bool executeCLI(const CLIOptions &commandOptions, int argc, char *argv[], const partition::PartitionOptions &partitionOptions,
		PackingOptions packingOptions, const SolutionMethodKPartition & solutionMethod) {
  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("num-partitions"	, partitionOptions.inputK						, UNSET_INT);
    optionRequired("file"						, partitionOptions.inputFilename		, UNSET_STRING);

    if (argc == 1 || vm.count("help")) {
      cout << helpOut.str();
      return true;
    }
  } catch(std::exception& e) {
    cout << helpOut.str()
         << "*** ERROR: " << e.what() << " ***"<< endl << endl;

    if (solutionMethod == RNP) {
    	cout << "*** NOTE: Make sure to run \"ulimit -s unlimited\" before running rnp! ***" << endl << endl;
    }
    return true;
  }
  return false;
}

// ================================================================
// Read from command line. Returns true if the program should
// exit due to help being asked for or some error, false otherwise.
// ================================================================
bool readCommandLine(int argc, char *argv[], partition::PartitionOptions &partitionOptions,
										 PackingOptions &packingOptions, SolutionMethodKPartition solutionMethod) {

	CLIOptions commandOptions(partitionOptions,packingOptions,solutionMethod);
	return executeCLI(commandOptions, argc,argv,partitionOptions, packingOptions, solutionMethod);
}



}
