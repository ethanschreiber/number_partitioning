/*
 * Utils.cpp
 *
 *  Created on: Mar 12, 2014
 *      Author: ethan
 */


#include "Utils.hpp"

// ************************
// String utility functions
// ************************

string getDirOfFile(string filename) {
  string fullpath(filename);
  size_t last = fullpath.find_last_of('/');
  return fullpath.substr(0,last+1);
}

string getFilename(string fullFilename) {
  string fullpath(fullFilename);
  size_t last = fullpath.find_last_of('/');
  return fullpath.substr(last+1,fullpath.size());
}

string addTrailingSlash(string dirname) {
  if (dirname[(dirname.size() - 1)] != '/') {
    dirname.append("/");
  }
  return dirname;
}

string getLastDir(string fullFilename) {
  string dir = getDirOfFile(fullFilename);
  dir = dir.substr(0,dir.size()-1);
  return getFilename(dir);
}

string removeTrailingSlash(string dirname) {
  if (dirname[(dirname.size() - 1)] == '/') {
    dirname.resize(dirname.size()-1);
  }
  return dirname;
}


bool startsWith(string s,string prefix) {

  // Note: strtncmp returns 0 if the first prefix.size() elements
  // are the same, hence the !

  return ((s.size() >= prefix.size()) &&
          (!strncmp(s.c_str(),prefix.c_str(),prefix.size())));
}

bool endsWith(string s,string suffix) {

  // Note: strtncmp returns 0 if the first prefix.size() elements
  // are the same, hence the !

  return ((s.size() >= suffix.size()) &&
          (0 == s.compare(s.length() - suffix.length(), suffix.length(),suffix)));
}

// **************************
// Filename utility functions
// **************************

bool isDir(const string &inputFilename) {
  struct stat s;
  if( stat(inputFilename.c_str(),&s) == 0 ) {     // If we can open file
    return s.st_mode & S_IFDIR;
  } else {
      cout << "Error opening [" << inputFilename << "]" << endl;
      exit(0);
  }
  return false;
}

// TODO: THis really belongs at the top level, it is for either packing or partitioning
// Helper function, you probably want to change this if you need to add options
static string getOutputFilenameHelper(const string filenameBase, const string suffix, const PackingOptions &packingOptions) {

  std::ostringstream outputFilename;
  outputFilename << filenameBase << suffix;

  // This is only relevant for binary-search bin packing  methods(BSBC or BELOV_BCP)
  if (packingOptions.useSchroeppelShamir) { outputFilename << SS_STRING; }
  if (packingOptions.classicSort)     	  { outputFilename << CLASSIC_SORT_STRING; }
  if (packingOptions.removePairs)         { outputFilename << PAIRS_STRING;  }
  if (packingOptions.useLDS)              { outputFilename << LDS_STRING; }

  // This is only relelvant if not using SS and not using classic Search
  if (!packingOptions.useSchroeppelShamir && (suffix.compare("_bsbc") == 0)) {
    outputFilename << "_IE" << packingOptions.bufferSize;
  }
  outputFilename << ".txt";
  return outputFilename.str();
}


string getOutputFilename(const string &inputFilename, const string suffix,const PackingOptions &packingOptions, int numPartitions /*= UNSET_INT */) {
  std::ostringstream outputFilename;

  if( isDir(inputFilename)) {
    DIR *dir = opendir (inputFilename.c_str());     // Open the directory

    if (dir != NULL) {                   // Make sure the directory can be opened

      if (numPartitions == UNSET_INT) {                             // For BinPacking
        outputFilename << addTrailingSlash(inputFilename)
                       << "output_all_";
      } else {
        outputFilename << addTrailingSlash(inputFilename)   // For NumberPartitioning
                       << "output_all_" << std::setfill('0') << std::setw(2) << numPartitions;
      }

    } else {
      cout << "Could not open directory [" << dir << "]" << endl;
      exit(0);
    }
  } else {
  	if (numPartitions == UNSET_INT) {
      outputFilename << getDirOfFile(inputFilename)
                     << "output_"
                     << getFilename(inputFilename);
  	} else {
			outputFilename << getDirOfFile(inputFilename)
										 << "output_"
										 << numPartitions << "_"
										 << getFilename(inputFilename);
  	}
  }

  return getOutputFilenameHelper(outputFilename.str(), suffix, packingOptions);
}



// skips "." and ".."  (current dir and parent dir)
// skips all strings beginning with "output_"
bool isInputFile(char *filename) {
  const char outputPrefix[8] = "output_";     // ignore files with this prefix
  return  (filename[0] != '.' &&                                     // Skip the current and parent directory
            strncmp(filename, outputPrefix, strlen(outputPrefix)));   // Skip the output files
}

vector<string> getInputFilenames(const string &inputFilename) {
  vector<string> filenames;   // The filenames to return


  string outputFilename;
  struct dirent *ent;                         // Dir Entry Struct
  if(isDir(inputFilename)) {                   // If the inputFilename is a directory
    DIR *dir = opendir (inputFilename.c_str());       // Open the directory
    while ((ent = readdir (dir)) != NULL) {  // read the files
      if (isInputFile(ent->d_name)) {         // If the file is a valid input file
        string filename = addTrailingSlash(inputFilename);
        filename.append(ent->d_name);
        filenames.push_back(filename);
      }
    }
  } else {
    string x(inputFilename);

    filenames.push_back(x);
  }

  std::sort(filenames.begin(),filenames.end());
  return filenames;
}

size_t countLines(const string filename, string prefix) {
  size_t numLines = 0;
  string line;
  std::ifstream inFile(filename);

  if (inFile) {
    while (std::getline(inFile, line)) {

      if (startsWith(line,prefix)) {
        ++numLines;
      }
    }
  }
  return numLines;
}
