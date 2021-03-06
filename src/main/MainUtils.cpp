/*
 * MainUtils.cpp
 *
 *  Created on: Mar 14, 2014
 *      Author: ethan
 */

#include "MainUtils.hpp"
#include <iostream>

using std::cout;
using std::endl;




// =====================================
// Read the problems from the input file
// =====================================
std::vector<partition::PartitionProblem> getProblems(const partition::PartitionOptions &partitionOptions) {
	std::vector<partition::PartitionProblem> problems;
	std::ifstream inFile(partitionOptions.inputFilename);
	std::string lineString;

	while (std::getline(inFile,lineString)) {
		std::stringstream lineStream(lineString);
		int idx;

		uint64_t element;
		std::vector<uint64_t> elements;
		lineStream >> idx;
		std::ostringstream problemName;
		problemName << "problem" << idx;
		uint64_t sum=0;
		while (lineStream >> element) {
			sum+=element;
			elements.push_back(element);
		}

		problems.push_back(partition::PartitionProblem(elements.size(),partitionOptions.inputK,&elements[0],
																					 sum, problemName.str()));
	}

	return problems;
}


std::vector<uint64_t> getSolutionValues(const string filename) {
	std::vector<uint64_t> solutionValues;
	std::ifstream inFile(filename);
	std::string lineString;

	while (std::getline(inFile,lineString)) {
		std::stringstream lineStream(lineString);
		string problemName;
		double time;
		uint64_t solutionValue;

		lineStream >> problemName >> time >> solutionValue;

		solutionValues.push_back(solutionValue);
	}
	return solutionValues;
}



// Sets buf and outFile
// Opens the output file, checks to see how many problems have been solved
// Either opens a new file or opens an existing one. Sets the streambuf to the proper point
// If only doing one file, prints to standard out instead
void openOutputFile(std::streambuf *& buf, std::ofstream &outFile, const string outputFilename, const int startIdx,
										const int problemIdx,	const size_t numProblems){

	{	// If in parallel, this prints the whole cout
		std::ostringstream outStream;
		outStream << "Output Filename: " << outputFilename << endl;
		cout << outStream.str();
	}

	if (problemIdx == UNSET_INT) {
		std::ostringstream os;
		os << "File: " << outputFilename << " ";

		if (startIdx >= (int) numProblems) {
			os << "[File already complete!]";
		} else if (startIdx > 0) {
			os << "[Appending from problem " << startIdx << ".]";
			outFile.open(outputFilename.c_str(), std::ios_base::app | std::ios_base::out);  // open in append
		} else {
			os << "[Starting from problem 0.]";
			outFile.open(outputFilename.c_str(), std::ios_base::trunc | std::ios_base::out); // truncate file first
		}
		buf = outFile.rdbuf();
		cout << os.str() << endl;
	} else {    // If just one file, use cout, no file necessary
		buf = std::cout.rdbuf();
	}

}
