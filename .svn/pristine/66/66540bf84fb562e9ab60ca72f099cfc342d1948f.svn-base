/*
 * convertTxtToBPA.cpp
 *
 *  Created on: Sep 18, 2012
 *      Author: ethan
 */


#include <iostream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdint.h>

using std::cout;
using std::endl;
using std::vector;
using std::string;


const string DIR_NAME("../belov_dat/");
const string FILENAME = DIR_NAME + "hard28";
int main(int argc, char *argv[]) {

  std::ifstream inFile(FILENAME.c_str());

  string name;
  int problemNum;
  int numEntries;
  int capacity;
  int entry,count;


  while (inFile) {
    string line;
    vector<int> numbers;
    int64_t sum = 0;
    char buf[255];
    std::getline(inFile, line);
    if(line.size() == 0) {
      break;
    }
    std::istringstream in(line);
    in >> buf;
    in >> problemNum;
    inFile >> numEntries;
    inFile >> capacity;

    cout << "Line       : [" << line << "]" << endl
         << "Problem Num: " << problemNum << endl
         << "Num Entries: " << numEntries << endl
         << "Capacity   : " <<  capacity << endl;

    for (int i=0;i<numEntries;i++) {
      inFile >> entry >> count;
      for (int j=0;j<count;j++) {
        numbers.push_back(entry);
        sum+= entry;
      }
    }
    std::ostringstream outputFilename;
    outputFilename << DIR_NAME << "problem" << problemNum << ".bpa";

    cout << "Filename: " << outputFilename.str() << endl;
    cout << "Num Entries: " << numbers.size() << endl << endl;
    inFile.ignore(255,'\n');

    std::ofstream outFile(outputFilename.str().c_str());

    outFile << "problem" << problemNum << endl
            << capacity << " " << numbers.size() << " " << sum << endl;
    for (size_t i=0;i<numbers.size();i++) {
      outFile << numbers[i] << endl;
    }
    outFile.close();
  }
  inFile.close();


}
