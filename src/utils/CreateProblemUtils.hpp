/*
 * CreateProblemUtils.hpp
 *
 *  Created on: Jan 18, 2013
 *      Author: ethan
 */

#ifndef CREATEPROBLEMUTILS_HPP_
#define CREATEPROBLEMUTILS_HPP_
#include <sys/stat.h>

void mkdir(string dir) {
  umask(0);               // Set the umask for mkdir to have correct permissions
  const char SEPARATOR[2] = "/";

  size_t separatorIdx = dir.find(SEPARATOR,0,1);
  while (separatorIdx!=string::npos) {

    mkdir(dir.substr(0,separatorIdx).c_str(), 0776);
    separatorIdx=dir.find(SEPARATOR,separatorIdx+1,1);
  }
}



#endif /* CREATEPROBLEMUTILS_HPP_ */
