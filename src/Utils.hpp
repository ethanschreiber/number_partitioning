/*
 * Utils.hpp
 *
 *  Created on: Mar 12, 2014
 *      Author: ethan
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "globals.hpp"
#include "pack/PackingUtils.hpp"
#include "partition/PartitionUtils.hpp"
#include <iostream>
#include <ctime>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <vector>


using std::cout;
using std::endl;
using std::string;
using std::vector;



// ============
// Timing Class
// ============
class SimpleTimer {
private :
  clock_t start;

public :
  SimpleTimer() : start(clock()) { }

  void reset() {
    start = clock();
  }
  double timeElapsed() {
    return (double(clock()) - double(start)) / CLOCKS_PER_SEC;
  }
};



// ************************
// String utility functions
// ************************

// Given a filename, removes everything after the last trailing slash
// and returns the result
string getDirOfFile(string filename);

// Given a full filename, removes the directory info and just
// return the filename
string getFilename(string fullFilename);

// Adds a trailing '/' to a dirname if it does not end with one
string addTrailingSlash(string dirname);

// /a/b/c/d/  <- Returns d
string getLastDir(string fullFilename);

// Removes a trailing '/' to a dirname if it ends with one
string removeTrailingSlash(string dirname);

// Returns true if s starts with prefix
// false otherwise
bool startsWith(string s,string prefix);

// Returns true if s eds with suffix
// false otherwise
bool endsWith(string s,string suffix);

template< typename T >
string join(const T v[], const size_t num, string delimiter=" ") {
  std::stringstream ss;
  for(size_t i = 0; i < num; ++i) {
    if(i != 0)
      ss << delimiter;
    ss << v[i];
  }
  return ss.str();
}

template< typename T >
string join(const vector<T> &v,string delimiter=" ") {
  return join(&v[0],v.size(),delimiter);
}


// **************************
// Filename utility functions
// **************************

bool isDir(const string &inputFilename);

// Include numPartitions for MultiwayNumberPartitioning, leave out for BinPacking
string getOutputFilename(const string &inputFilename, const string suffix, const PackingOptions &packingOptions, int numPartitions);

bool isInputFile(char *filename);
vector<string> getInputFilenames(const string &inputFilename);

// Count the number of lines in  filename starting with prefix
size_t countLines(const string filename, const string prefix);


#endif /* UTILS_HPP_ */
