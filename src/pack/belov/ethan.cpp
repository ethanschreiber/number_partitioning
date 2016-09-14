#include "../BinCompletionUtils.hpp"
#include "stdafx.h"
#include "solver.h"
#include "bcp.h"  // uses problem.h
#include "bbcuts.h"
#include "probl_csp1.h"
#include "probl_csp2.h"
#include "probl_cp22.h"
#include "probl_pmp1.h"
#include "bcp.h"
#include "bcp2.h"
#include "raster.h"
#include "lasthdr.h"


////////////////////////////////////
#include <ctype.h>

// --------
// main.cpp
// --------
//Ethan: Alternative to main()
void executeBelovBCP( const BinPackingProblem &problem,ProblemStats &stats) {
  bool usePerfectPairs = false;
  ss::Solver slv;
  ss::Solver::RWOptions(); // even if no files

  BinPackingProblem processedProblem(problem.capacity);
  PairVector perfectPairs;

  clock_t start = clock();
  if (usePerfectPairs) {
    stats.numPerfectPairs = removePairs(problem,processedProblem,perfectPairs);
    slv.processStruct(processedProblem,stats);
  } else {
    stats.numPerfectPairs = 0;
    slv.processStruct(problem,stats);
  }
  clock_t finish = clock();

  stats.numBins    += stats.numPerfectPairs;
  stats.lowerbound += stats.numPerfectPairs;
  stats.sum = problem.sum;
  stats.time = (double(finish)-double(start))/CLOCKS_PER_SEC;

}


using namespace COMMON_NAMESPACE__;
SS_BEGIN_NAMESPACE__
// Ethan ProblemLoader for BinPackingProblem struct
int ProblemLoader::Read(const BinPackingProblem &problem) {
  infile = "nofile.dat";
  strcpy(prNm, problem.problemName.c_str());
  inst = 1;
  pr = new CSP1(penv, infile,inst,prNm);
  return pr->Read(problem);
}

// Ethan: Alternative to processFile
void Solver::processStruct(const BinPackingProblem &problem,ProblemStats &stats) {
  Env env;

  ProblemLoader pl(env,NULL,1);
  int readRes   = pl.Read(problem);

  auto_ptr<BCP> bcp(pl.GetSolver());
  bcp->firstinstance = true;

  bcp->Run();


  //cout << "# Bins: " << bcp->gub << endl;
  stats.numBins = bcp->gub;
  stats.lowerbound = bcp->glb0;
  stats.numNodes = bcp->iterAll;

}

// ---------------
// probl_cpp22.cpp
// ---------------
// Ethan: Not used but implemented to avoid pure virtual error
int CP22::Read(const BinPackingProblem &problem) {
  cout << endl << "ERROR: CP22::READ UNIMPLEMENTED, EXITING!" << endl << endl;
  exit(0);
}

// ---------------
// probl_csp2.cpp
// ---------------
// Ethan: Not used but implemented to avoid pure virtual error
int CSP2::Read(const BinPackingProblem &problem) {
  cout << endl << "ERROR: CSP2::READ UNIMPLEMENTED, EXITING!" << endl << endl;
  exit(0);
}

// ---------------
// probl_csp1.cpp
// ---------------

// Ethan: For reading from BinPackingProblem struct
int CSP1::Read(const BinPackingProblem &problem) {
  L0 = problem.capacity;
  m = problem.N;      // TODO: This is not quite right if there are duplicates
  m0 = problem.N;

  pc0.resize(m0);

  for (int i=0;i<m0;i++) {
    pc0[i].l = problem.S[i];    // The item size
    pc0[i].b = 1;               // # of copies of this item, TODO: Fix this, there might be duplicates
  }

  InitProblem();
  return 0;
}
SS_END_NAMESPACE__
