// solver.h: Schnittstelle f�r die Klasse Solver.
//
//////////////////////////////////////////////////////////////////////


#include "bcp.h"
#include "../PackingUtils.hpp"    // Ethan: So we can use BinPackingProblem struct

SS_BEGIN_NAMESPACE__

class ProblemLoader {
public:
  Env * penv;
  const char * infile;
  int inst;
  char prNm[100];
  Problem * pr;

  ProblemLoader(Env & e,const char *in,const int i)
    : penv(&e), infile(in), inst(i) { }

  int Read(istream &);
  // Ethan: Added function to read from struct
  int Read(const BinPackingProblem &problem);

  BCP * GetSolver();
};//____________________________________________________

class Solver : public MyApp
{
public:
	virtual void CalculateTimeUnit();
	virtual void ProcessFile(const char * fln);

	// Ethan: Alternative to ProcessFile
	  virtual void processStruct(const BinPackingProblem &problem,ProblemStats &stats);
  virtual void ProcessPr6(ProblemLoader *);
  virtual void ProcessPr7(ProblemLoader *);
  static void RWOptions();
  static void WriteOptions();
	//Solver();
  virtual ~Solver() { }
  static int problemType;
  static int instFirst, instLast;
  static int fTestParams;

  static int nTries, nRounds, nParts;
  static double splitTimeLimit;
  static int minL, maxL;
  Random rndGrp; // for grouping items
  static double minZ; // to access from PMP
protected:
//	virtual void DoneFile()=0;
//	virtual void InitFile()=0;

  static opt::OptSection opt;
  static opt::OptContainer Options();
  static double outputLevel;
};

SS_END_NAMESPACE__

