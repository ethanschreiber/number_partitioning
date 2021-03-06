/*
 * SubsetSumMain.hpp
 *
 *  Created on: Oct 1, 2012
 *      Author: ethan
 */

#ifndef SUBSETSUMMAIN_HPP_
#define SUBSETSUMMAIN_HPP_
#include "Utils.hpp"
#include "pack/PackingUtils.hpp"
#include "ss/SubsetSum.hpp"
#include "ss/SchroeppelShamir.hpp"
#include "ss/extended/InclusionExclusion.hpp"
#include "ss/extended/Extended_Schroeppel_Shamir.hpp"
#include "ss/extended/Extended_Horowitz_Sahni.hpp"


#include <algorithm>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>   // For ceil
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <vector>


using std::vector;
using std::string;
using std::ofstream;

using std::cout;
using std::endl;


string getOutputFilename(const string &inputFilename,
                         const double lb, const double ub) {
  std::ostringstream outputFilename;

  if( isDir(inputFilename)) {
    DIR *dir = opendir (inputFilename.c_str());     // Open the directory

    if (dir != NULL) {                   // Make sure the directory can be opened

      outputFilename << addTrailingSlash(inputFilename)
                     << "output_all_subset_sum_";
    } else {
      cout << "Could not open directory [" << dir << "]" << endl;
      exit(0);
    }
  } else {
    outputFilename << getDirOfFile(inputFilename)
                   << "output_subset_sum_"
                   << getFilename(inputFilename);
  }

  outputFilename << lb << "_" << ub << ".txt";

  return outputFilename.str();
}

#endif /* SUBSETSUMMAIN_HPP_ */
