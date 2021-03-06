/*
 * ProcessExperimentsUtils.hpp
 *
 *  Created on: Jan 24, 2013
 *      Author: ethan
 */

#ifndef PROCESSEXPERIMENTSUTILS_HPP_
#define PROCESSEXPERIMENTSUTILS_HPP_

#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include "../pack/PackingUtils.hpp"
#include "../Utils.hpp"
// ****************************************************************************
// Given the ./dat/ directory, returns all experimental directories within
// that directory in a sorted vector. i.e., all dirs ending with DIR_SUFFIX
// ****************************************************************************
vector<string> getSubDirs(const char *datDirName) {
  vector<string> subdirNames;

  // Read the experiment directories, fill in experimentDirNames vector
  DIR *datDir = opendir (datDirName);         // Open the directory
  struct dirent *ent;                         // Dir Entry Struct

  while ((ent = readdir (datDir)) != NULL) {  // read the files

    string dirName(ent->d_name);

    if (!endsWith(dirName,".")) {
      string subdirName = addTrailingSlash(datDirName);
      subdirName.append(addTrailingSlash(dirName));
      subdirNames.push_back(subdirName);
    }
  }

  // Sort the directories alphabetically
  std::sort(subdirNames.begin(),subdirNames.end());

  return subdirNames;
}


#endif /* PROCESSEXPERIMENTSUTILS_HPP_ */
