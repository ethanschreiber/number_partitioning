/*
 * ProcessSolution.cpp
 *
 *  Created on: Sep 27, 2012
 *      Author: ethan
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdint.h>
using std::cout;
using std::endl;
using std::string;
using std::vector;
int main(int argc, char *argv[]) {
  std::cout.imbue(std::locale("")); // For printing 1,000 instead of 1000

  if (argc != 2) {
   cout << endl
        << " Usage: " << argv[0] << " [input file] " << endl << endl;
   exit(0);
  }
  char *& inputFilename = argv[1];

  std::ifstream inFile(inputFilename);

  int64_t N,K,totalLB,totalUB,isOptimal;

  inFile >> K >> N >> totalLB >> totalUB >> isOptimal;

  cout << "N         : " << N << endl
       << "K         : " << K << endl
       << "LB        : " << totalLB << endl
       << "UB        : " << totalLB << endl
       << "Is Optimal: " << (isOptimal ? "yes" : "no") << endl << endl;



  inFile.ignore(255,'\n');

  vector<vector<int> > indicesMatrix;
  vector<vector<int> > numbersMatrix;
  vector<int> sums;
  for (int i=0;i<K;i++) {
    string line1;
    string line2;
    int sum;

    vector<int> v1;
    vector<int> v2;
    std::getline(inFile,line1);
    std::getline(inFile,line2);
    inFile >> sum;
    inFile.ignore(255,'\n');

    std::istringstream iss1(line1);
    std::istringstream iss2(line2);

    copy( std::istream_iterator <int> ( iss1 ),
          std::istream_iterator <int> (),
          std::back_inserter( v1 ) );
    indicesMatrix.push_back(v1);


    copy( std::istream_iterator <int> ( iss2 ),
          std::istream_iterator <int> (),
          std::back_inserter( v2) );
    numbersMatrix.push_back(v2);

    sums.push_back(sum);

  }

  //  for (size_t i=0;i<numbersMatrix.size();i++) {
  //    cout << std::setw(2) << i << ": (Sum = " << sums[i] << ") ";
  //    int sum =0;
  //    for (size_t j=0;j<numbersMatrix[i].size();j++) {
  //      cout << numbersMatrix[i][j] << " ";
  //      sum += numbersMatrix[i][j];
  //    }
  //
  //    if (sum != sums[i]) {
  //      cout << "SUMS NOT EQUAL: " << sum << " != " << sums[i];
  //    }
  //
  //    cout << endl;
  //  }

    for (size_t count=0;count<numbersMatrix.size();count++) {
      cout << std::setw(2) << count << ": (Sum = " << sums[count] << ") ";
      int sum =0;

      size_t maxRow = 0;
      int64_t maxValue = numbersMatrix[0][0];
      for (size_t row=1;row<numbersMatrix.size();row++) {
    	  if (numbersMatrix[row][0] > maxValue) {
    		  maxValue = numbersMatrix[row][0];
    		  maxRow = row;
    	  }
      }

      for (size_t j=0;j<numbersMatrix[maxRow].size();j++) {
        cout << numbersMatrix[maxRow][j] << " ";
        sum += numbersMatrix[maxRow][j];
      }
      numbersMatrix[maxRow][0] = -1;

      if (sum != sums[maxRow]) {
        cout << "SUMS NOT EQUAL: " << sum << " != " << sums[count];
      }

      cout << endl;
    }


}
