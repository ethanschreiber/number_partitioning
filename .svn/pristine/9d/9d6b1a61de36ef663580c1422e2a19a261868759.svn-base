#ifndef __BB_H
#define __BB_H

// How data storage ? 0-based arrays. Auto-dimensioning
// 1. Level cut 2. Dual simpl 3. +Cuts. Constraint pool

#include "problem.h"

SS_BEGIN_NAMESPACE__

/* Before using the class, call Reallocate(m,mc).
 Before Run(), input vars must be set.
  */
class BB {
public:
  struct Piece {
    size l;
    int
      b,
      i0;       // the original index before sorting
    double
      d,        // the simplex mult. or weight
      q,
//      dsa,
      qma,      // q max after (inclusive)
      w;        // weight for sorting
    size
      lma      // l min after (not inclusive)
//     , lsa
      ;      // l sum after
    bool operator > (const Piece & pc2) const
    { return w > pc2.w; } // for sorting (w!!)
  };//__________________________________________________
  typedef Vector<Piece> PieceContainer;
  //////////////////////////////////////////////////////
  /// INPUT:: (not to be modified inside)
  //____________________________________________________
/// Number of piece types //, addi constr:
  int m;
  size L; // Stock length
  PieceContainer pieces__; // Each knows l,b,d
  Pool<Column> * forbidden__;
  bool fHeur; // whether called from heuristic, then no forbidden
/// The additional constraints part of the column, fixed:
  //valarray<addi_float> addi__;
  //valarray<double> dAddi__; // their weights
/// Initial, eg incumbent, lower bound:
  double zLowerInitial;
/// A priori lower bound, eg for primal simplex or
/// obtained from Lagrange relaxation.
/// E.g. CSP1: 1
  double zMin;
/// The tolerance by which a better solution must be >=:
  double eps; // =smth*bb_eps
//  static bool fRandomizeWeights;
// + How -- params
  static double nStepsTooMuch;
  static bool fEvenIfNotFound;
  static double outputLevel;
  //////////////////////////////////////////////////////
  /// OUTPUT::
  //____________________________________________________
/// Result flag. =true e.g. if (z > zMin).
  bool found;
  bool fETerm;
/// The value of the best found solution
  double z;
/// Solutions satisf. the lower bound:
  ColSet *colsetRes;
/// The best sol.
  Column *colBest; // gives also initial info in 2D
  //////////////////////////////////////////////////////
  /// STATISTICS::
//  mystat::Accumulator
    // Init: some for 1 problem, some for 1 iteration,
    // some for all problems!
  //////////////////////////////////////////////////////
private:
  //double zfConstTerm;
  static opt::OptSection opt;
  //////////////////////////////////////////////////////
  /// Methods::
  //____________________________________________________
public:
  BB() :eps(1e-8), colsetRes(NULL), colBest(0),
    m(0), L(0), z(0), zLowerInitial(0), zMin(0), fHeur(0), forbidden__(NULL) { }
  virtual ~BB() { }
  virtual void Reallocate(int m);
// Precondition: All input vars set!
  virtual void Run(); // need this ??
  // Randomization? Restarts? Time control ?
  // Perturbation: slightly change weights
protected:
  virtual void Optimize();
  void SignalBetterSolution();
  void RandomizePieceWeightRatio(Piece *pc);
  static opt::OptContainer Options();
  virtual void Print();
  void DoStatistics(int what, ostream &);
}; // BBCuts
  //____________________________________________________



SS_END_NAMESPACE__

#endif // __BB_H
