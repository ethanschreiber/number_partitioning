/*
 * ProcessExperiments.cpp
 *
 *  Created on: Aug 16, 2012
 *      Author: ethan
 */

#include "../pack/PackingUtils.hpp"
#include "ProgramOptionsUtils.hpp"
#include "../globals.hpp"
#include "ProcessExperimentsUtils.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cmath>

using namespace boost;
namespace po = boost::program_options;
using std::cout;
using std::endl;
using std::vector;
using std::map;


const int DECIMAL_PRECISION = 4;
const int NUM_FILES = 25;    // TODO: Read this in

const int NUM_P = 3;    // Number of different precisions (orders of magnitude}
const int P_VALUES[NUM_P] = {4,6,15};

//const int NUM_P = 2;    // Number of different precisions (orders of magnitude}
//const int P_VALUES[NUM_P] = {4,15};

//// This is for partitioning 45
//const int NUM_METHODS = 2;  // Number of different methods
//const string METHOD_NAMES_FULL[NUM_METHODS] = {"BSIBC", "BSBCP"};
//const string METHOD_FILE_NAMES[NUM_METHODS] = {"output_all__bc_IE50.txt",
//                                               "output_all__belov_bcp.txt"};


// This is for packing 100
//const int NUM_METHODS = 3;  // Number of different methods
//const string METHOD_NAMES_FULL[NUM_METHODS] = {"Improved BC", "Original BC", "BCP"};
//const string METHOD_FILE_NAMES[NUM_METHODS] = {"output_all__bc_IE50.txt",
//                                               "output_all__bc_oldsort_IE1000000.txt",
//                                               "output_all__belov_bcp.txt"};

//
//// This is for packing 45
const int NUM_METHODS = 2;  // Number of different methods
const string METHOD_NAMES_FULL[NUM_METHODS] = {"IBC",  "BCP"};
const string METHOD_FILE_NAMES[NUM_METHODS] = {"output_all__bc_IE50.txt",
                                               "output_all__belov_bcp.txt"};


// This is for partitioning 45 with only 1 improvement
//const int NUM_METHODS = 2;  // Number of different methods
//const string METHOD_NAMES_FULL[NUM_METHODS] = {"Sort",  "Dynamic"};
//const string METHOD_FILE_NAMES[NUM_METHODS] = {"output_all__bc_IE1000000.txt",
//                                               "output_all__bc_oldsort_IE50.txt"};

//// This is for packing 100 with only 1 improvement
//const int NUM_METHODS = 2;  // Number of different methods
//const string METHOD_NAMES_FULL[NUM_METHODS] = {"Sort",  "Dynamic"};
//const string METHOD_FILE_NAMES[NUM_METHODS] = {"output_all__bc_IE1000000.txt",
//                                               "output_all__bc_oldsort_IE50.txt"};

const int COL_WIDTH=10;

// ****************************************************************************
// Main
// ****************************************************************************
int main(int argc, char *argv[]) {

  int N = UNSET_INT;

  std::ostringstream OUTPUT_FILENAME_SS;
  std::ostringstream DAT_DIR_SS;

  bool isPartitioning = false;
  po::options_description commandOptions("Options");
  commandOptions.add_options()
  ("help,h"        , "produce help message")
  ("partitioning,p", po::value(&isPartitioning)->zero_tokens(),"Use this for partitioning as opposed to packing")
  ("numElements,n"           , po::value< int >(&N), "(Default=100) The number of input items per experiment");

  std::ostringstream helpOut;
  helpOut << endl << "Usage: " << argv[0] << " [options]" << endl << endl << commandOptions << endl;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, commandOptions), vm);
    po::notify(vm);

    optionRequired("numElements",N,UNSET_INT);
    if (argc == 1 || vm.count("help")) { cout << helpOut.str();  return 0;  }
  } catch(std::exception& e) {
    cout << helpOut.str() << "*** ERROR: " << e.what() << " ***"<< endl << endl;
    return 1;
  }

  if (isPartitioning) {
    OUTPUT_FILENAME_SS << "/home/ethan/workspace/ijcai13/partitioning_table_" << N << ".tex";
    DAT_DIR_SS         << "/home/ethan/workspace/BinPacking/dat/partitioning/" << N << "/";
  } else {
    OUTPUT_FILENAME_SS << "/home/ethan/workspace/ijcai13/packing_table_" << N << ".tex";
    DAT_DIR_SS         << "/home/ethan/workspace/BinPacking/dat/packing/" << N << "/";
  }

  const string OUTPUT_FILENAME = OUTPUT_FILENAME_SS.str();
  const string DAT_DIR         = DAT_DIR_SS.str();

  std::ostringstream latex;

  // ============
  // Latex header
  // ============
  latex << "\\begin{tabular}{|rll |";

for (int mIdx = 0;mIdx < NUM_METHODS;mIdx++) {
  latex << "|";
    for (int pIdx=0;pIdx<NUM_P;pIdx++) {
      latex << " r ";
    }
}
  latex << "|}" << endl;
  string header = (isPartitioning) ? "Magnitude (M)" : "Bin Capacity (C)";
  latex << "\\cline{4-" << ((NUM_METHODS * NUM_P) + 3) << "} \\multicolumn{3}{c|}{N = \\textbf{" << N << "}} & "
        <<  "\\multicolumn{" << NUM_METHODS * NUM_P << "}{|c|}{\\textbf{" << header << "}} \\\\" << endl;
  latex << "\\hline" << endl;
  latex << "& &   ";
  for (int i=0;i<NUM_METHODS;i++) {
    latex << "& \\multicolumn{" << NUM_P << "}{c|}{\\textbf{" <<  METHOD_NAMES_FULL[i] << "}}";
  }
  latex << "\\\\" << endl;

  if (isPartitioning) {
    latex << "\\multicolumn{2}{|c}{\\textbf{Sample Range}} & \\textbf{k}";
  } else {
    latex << "\\multicolumn{2}{|c}{\\textbf{Sample Range}} & \\textbf{E}$[B]$ ";
  }
  for (int mIdx = 0;mIdx < NUM_METHODS;mIdx++) {
    for (int pIdx=0;pIdx<NUM_P;pIdx++) {
      std::ostringstream pOut;
      pOut << " $10^{" << P_VALUES[pIdx] << "}$";
      latex << "&" << std::right << std::setw(COL_WIDTH) << pOut.str() << " ";
    }
  }
  latex << "\\\\" << endl
        << "\\hline" << endl;


  // read in fValues


  vector<string> fValues;
  vector<int> expectedValues;
  {
    const double N_2 = ((double) N / (double) 2);
    std::ostringstream fValueFilename;
    fValueFilename << DAT_DIR << "10_" << P_VALUES[0] << "/";

    fValues = getSubDirs(fValueFilename.str().c_str());
    for (size_t i=0;i<fValues.size();i++) {
      fValues[i] = getLastDir(fValues[i]);

      double F = atof(fValues[i].c_str());
      expectedValues.push_back((int) std::ceil(N_2 * F));
    }
  }

  // ===========================
  // Read from experimental dirs
  // ===========================

  for (size_t fIdx = 0;fIdx < fValues.size();fIdx++) {
    for (int methodIdx = 0;methodIdx<NUM_METHODS;methodIdx++) {

      if (methodIdx == 0) {

        string constant = isPartitioning ? "M" : "C";
        latex << std::left << std::setw(5) << "$[1,$&$" << fValues[fIdx] << "C]$&" << std::setw(5)
              << expectedValues[fIdx];
      } else {
        latex << std::setw(11) << "";
      }
      for (int pIdx = 0;pIdx < NUM_P;pIdx++) {
        std::ostringstream filename;
        filename << DAT_DIR << "10_" << P_VALUES[pIdx] << "/" << fValues[fIdx] << "/" << METHOD_FILE_NAMES[methodIdx];

        std::ifstream inFile(filename.str());

        double totalTime = 0;
        if (inFile) {   // If the file exists

          int count = 0;
          string line;

          while (std::getline(inFile, line)) {

            std::stringstream  lineStream(line);
            string        cell;
            string problemName;
            string timeString;
            double time;

            if (isPartitioning) {
              std::getline(lineStream,problemName,' ');
              std::getline(lineStream,timeString,' ');
              time = ::atof(timeString.c_str());

            } else {
              std::getline(lineStream,problemName,',');
              std::getline(lineStream,timeString,',');
              time = ::atof(timeString.c_str());
            }
            totalTime += time;
            count++;
          }

          if (count == NUM_FILES) {
            totalTime = totalTime / (double) count;
          } else {
            totalTime = UNSET_DOUBLE;

          }
        } else {
          totalTime = UNSET_DOUBLE;
          cout << filename.str() << " does not exist" << endl;
        }


        latex << "&" << std::left << std::setw(COL_WIDTH);
        if (totalTime == UNSET_DOUBLE) {
          latex << std::fixed << std::setprecision(DECIMAL_PRECISION) << "-" << " ";
        } else {
          latex << std::fixed << std::setprecision(DECIMAL_PRECISION) << totalTime << " ";
        }

      }
      if (methodIdx == NUM_METHODS - 1) {
        latex << "\\\\";
      }
        latex << endl;
//      if (methodIdx == NUM_METHODS - 1) {
//            latex << "\\hline" << endl;
//      }

    }
  }

  // ============
  // Latex Footer
  // ============
  latex << "\\hline" << endl
        << "\\end{tabular}" << endl;


  cout << "Writing tex to " << OUTPUT_FILENAME << endl;
  std::ofstream outFile(OUTPUT_FILENAME);
    outFile << latex.str();
    outFile.close();

}
