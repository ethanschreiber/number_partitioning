/*
 * CLIW.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: ethan
 */

#ifndef CLI_HPP_
#define CLI_HPP_

#include "../partition/PartitionUtils.hpp"
#include "../utils/ProgramOptionsUtils.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
using std::string;


namespace po = boost::program_options;

// Specification
namespace cli {

/*************************************************
 * The options to pass to boost::program_options *
 *************************************************/
class CLIOptions : public po::options_description {
public :
	void addPartitionOptions(partition::PartitionOptions &partitionOptions) {
		add_options()
		  ("help,h"           , "produce help message")
		  ("num-partitions,k" , po::value< int >    (&partitionOptions.inputK)            , "The number of partitions. (required)")
		  ("file,f"           , po::value< string > (&partitionOptions.inputFilename)     , "(required) The input filename.")
		  ("problem-index,x"  , po::value< int >    (&partitionOptions.problemIdx)        , "The problem index from the input filename. If this is left out, all problems in the file are solved");
	}

	CLIOptions(partition::PartitionOptions &partitionOptions, PackingOptions &packingOptions,
					SolutionMethodKPartition solutionMethod) : po::options_description("Options") {

		addPartitionOptions(partitionOptions);


		if (solutionMethod == CIW || solutionMethod == CIW_LC) {
		  	add_options()
		      ("lc,l"             , po::value(&partitionOptions.isLowCardinality)             ->zero_tokens(),"Should we use the low cardinality version? (flag)")
		  		("num-cached-sets,n", po::value< size_t > (&partitionOptions.numSets)           ,"(DEFAULT=20000 for normal, 1000 for low cardinality) The number of cached sets to generate for each schroeppel and shamir run.");
		} else if (solutionMethod == MOFFITT || solutionMethod == RNP || solutionMethod == SNP || solutionMethod == RNP_2009 || solutionMethod == CGA_MW){
			// Nothing to add for MOF or RNP or SNP
		} else if (solutionMethod == BSBC || solutionMethod == BSBCP) {
		  add_options()
		  ("lower-bound,l"    , po::value< uint64_t >(&partitionOptions.minCapacity)   , "The lower bound for the binary search.")
		  ("upper-bound,u"    , po::value< uint64_t >(&partitionOptions.maxCapacity)   , "The upper bound for the binary search.")
		  ("buffer-size,b"    , po::value< int >(&packingOptions.bufferSize)            , "(Default=50) The size of the buffer for inclusion/exclusion search")
		  ("lds,d"            , po::value(&packingOptions.useLDS)                          ->zero_tokens(),"Should we use limited discrepency search (flag, include or don't include.")
		  ("ss,s"             , po::value(&packingOptions.useSchroeppelShamir)             ->zero_tokens(),"Should we use Schroeppel and Shamir to generate completions? (flag)")
		  ("ssbs,t"           , po::value(&packingOptions.useSchroeppelShamirBinarySearch) ->zero_tokens(),"Should we use Schroeppel and Shamir to generate possible values for the binary search? (flag)")
		  ("printsol,p"       , po::value(&packingOptions.printSolution)                   ->zero_tokens(),"Should we print the solution? (flag)")
		  ("verbose,v"        , po::value(&packingOptions.isVerbose)                       ->zero_tokens(),"Should we display verbose output? (flag)")
		  ("classic-sort,c"   , po::value(&packingOptions.classicSort)->zero_tokens(),"Should we use the classic korf bin completion sort?") ;

		} else {
			cout << "In CLI.hpp: Invalid solution method: " << solutionMethod << endl;
			exit(0);
		}

	}
};

bool readCommandLine(int argc,char *argv[],partition::PartitionOptions &partitionOptions,PackingOptions &packingOptions,
										SolutionMethodKPartition solutionMethod);
}



#endif /* CLI_HPP_ */
