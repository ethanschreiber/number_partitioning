// FILE: bcp_lp.cpp, solution of LPs in branch&cut&price
// Author: Gleb <Belov@math.tu-dresden.de>

#include "stdafx.h"
#include "bcp.h"
#include "solver.h"
#include "lasthdr.h"

/// !!! Unlike ABACUS, we don't generate cuts
/// if pricing time; more, we delete cuts

SS_BEGIN_NAMESPACE__

// Distinguish LP status

//bool BCP::StopLP()
//{ return TimeLimit() or LP::error == lp->status
//  /*or LP::unbounded == lp->status*/; }

void BCP::SetCurrentLPV(double v) // during LP
{ clb = pr->raster_ceil(clv=v); // for a restricted master
pr->SetCurrentLPBnd(clb);
pr->SetCurrentLPValue(clv); }
void BCP::SetLocalLPV(double v) { // in locla LP opt
  llb = pr->raster_ceil(llv=v);
  pr->SetLocalLPBnd(llb); // BUT CARE: when raster changes?
}
void BCP::UpdateLocalLPV() { // set from current (clv)
  llb = clb; llv = clv;
  pr->SetLocalLPBnd(llb); // BUT CARE: when raster changes?
}
void BCP::UpdateGlobalLPV() { // better lpv value
  glb = FMin(glb_others,llb); glv = FMin(glv_others,clv);
  if (glb!=pr->GetLPBnd())
    timeLPBest = timer.userTime();
  pr->SetLPBnd(glb); // BUT CARE: when raster changes?
    // when raster_ceil called with gub<inf, glb>-inf;
}
double BCP::GetLocalLagrBnd() {
//  if (lLv != pr->GetLocalLagrValue()) {
    lLv = pr->GetLocalLagrValue();
    lLb = pr->raster_ceil(lLv);
//  }
  return lLb;
}

int BCP::GetColIndex(int i) {
  if (cols[i].slackCut)
  { dbg_outn_(1,'x'); return -1; }
  return cols[i].j;
}
Column * BCP::GetCol(int i) {
  assert(i>=0 and i<cols.size());
  assert(not cols[i].slackCut);
  return GetMainCol(cols[i].j);
}
Column * BCP::GetMainCol(int i) {
  assert(0<=i and i<colpool.size());
  assert(colpool[i]);
  return colpool[i];
}
bool BCP::IsSlackCol(int i) // the cols of LP
{return cols[i].slackCut ? true: GetCol(i)->IsSlack();}


void BCP::InitLP() {
  pr->InitLP();
  fLPOpt = false;
  clb=clv=1e100; // current lp bnd/val
//  llbOld = llb; // ??
  dbg_outn_(2," LP:");
}

bool BCP::SolvePrimal() {
  fLPIntOpt = false;
  lp->SolvePrimal();
  lp->GetMultipliers(lpd); // even if infeas
  GetLPSolution();
  dbg_outn_(3," s");
  if (LP::opt != lp->status)
  { dbg_outn_(3," lpNO"); return false; }
  if (fInfeasBasic)
  { dbg_outn_(3," Inf"); return true; }
  SetCurrentLPV(GetLPValue());
  dbg_outn_(3," LP(opt)val: "<<clv);
  Rounding(2);
  return true;
}
bool BCP::SolveDual() {
  fLPIntOpt = false;
  lp->SolveDual();
  lp->GetMultipliers(lpd); // even if unbounded
  GetLPSolution();
  dbg_outn_(3," S");
  if (LP::opt != lp->status)
  { dbg_outn_(3," lpNO"); return false; }
  if (fInfeasBasic)
  { dbg_outn_(3," Inf"); return true; }
  SetCurrentLPV(GetLPValue());
  Rounding(3);
  return true;
}
bool BCP::SolveDual_GetMultsOnly() { // for pscosts
  lp->SolveDual();
  lp->GetMultipliers(lpd); // even if unbounded
  dbg_outn_(3," Sd");
  if (LP::opt != lp->status)
  { dbg_outn_(3," lpdNO"); return false; }
  return true;
}
void BCP::GetLPSolution() {
  int i;
  lp->GetSolution(&cstat,&lpx);
  if (LP::opt != lp->status) return;
  fInfeasBasic = 0;
  log_n_(3," LPv=" <<setprecision(16) << GetLPValue()
    << " fracpart="<<smfrac(GetLPValue())<<setprecision(6));
  for (i=0;i<lpx.size();++i) {
    if (lpx[i]<-1e-6) {
      lp->WriteModel("neg.lp");
      GetEnv().GetLog2() << "LP solution with neg vars:\n";
      PrintContSol();
      assertm(0,"\nNegative vars!!! See neg.lp && [input]__.txt");
    }
    if (lpx[i] > pr->GetXEps()) {
     bool fInf = false;
     if (not cols[i].slackCut)
      if (GetCol(i)->fInfeasible()) fInf = 1;
     if (-2 == cols[i].j) // infeas cut slack
      fInf = 1;
     if (fInf) {
      dbg_outn(3, "Infeas column "<<i<<" has value "<<lpx[i]);
      fInfeasBasic = 1;
     }
    }
    if (OUTP_LEV__ >= 3.5)
      if (lpx[i] && not cols[i].slackCut)
        log__(" x"<<GetColIndex(i)<<'='<<lpx[i]);
  }
}
// returns 0 if LP is infeas even with addi slacks
// -- OLD -- now: only updates coef
bool BCP::StartPhaseI() {
/*
  int i;
  for (i=0;i<cols.size();++i) {
    if (not cols[i].slackCut)
    if (GetCol(i)->fInfeasible()) {
//      lp->ChangeVarBnds(i,0,1e100);
      lp->SetObjCoef(i, maxObjCoef * M__);
    }
//    else
  //    lp->SetObjCoef(i,0); // for non-phaseI vars
  } */ // not changed
  log_n_(3," SPhaseI: infeas cols obj coef (NOT SET)="<<maxObjCoef * M__);
  if (OUTP_LEV__ >= 7) lp->WriteModel("PhaseI.lp");
/*  if (!SolveDual()) {
    EndPhaseI();
    return 0;
  }
  log_n_(3," SPhaseI started.");*/
  fInfeasColsIn = 1; // why here?
  return 1;
}
void BCP::EndPhaseI() {
  return; /*
  int i;
  for (i=0;i<cols.size();++i) {
    if (GetCol(i)->fInfeasible())
      lp->ChangeVarBnds(i,0,0);
    //else
      //lp->SetObjCoef(i,GetCol(i)->GetObj()); // for non-phaseI vars
  }
  fInfeasColsIn = 0;
  log_n_(3," EPhaseI");
//  if (!SolvePrimal()) { // no resolving.
*/
}
bool BCP::Price() { // Subgradient Optimization when no cuts
  dbg_outn_(3," _price ");
  int k=1;  // :: if (pr->DualBndEqualToLP())
  if (clb == llbLast)
  if (fSkipColGenAtLPBound) { // or value
    // fLPOpt = ? fabs(clv-llvLast)<1e-6;
    dbg_outn_(3," skipCG_");
    goto RetFalse;
  }
// PriceAgain:
  if (!GenCol(lpd)) { // RMP opt
    if (fInfeasBasic)
/*      if (maxObjCoef < 1e20// && M1__<M__max
      ) {
      maxObjCoef *= 10;
      log_n_(1," Inf. cols basic. max objCoef="<<maxObjCoef);
      for (int i=0;i<cols.size();++i)
      if (GetCol(i)->fInfeasible()) {
//      lp->ChangeVarBnds(i,0,1e100);
        lp->SetObjCoef(i, maxObjCoef * M__);
      }
      log_n_(3," SPhaseI: infeas cols obj coef="<<maxObjCoef * M__);
      if (!SolvePrimal())
        goto RetFalse;
      goto PriceAgain;
      // NOT ENDING PHASE I
      }
      else */ {
        log_n_(3," IX");
        lp->status = LP::infeas;
        goto RetFalse;
      }
/*    else
      if (fInfeasColsIn) {
        EndPhaseI();
        return true; // !!!
      }*/
    log_n_(3,'X');
    fLPOpt = true;
    goto RetFalse;
  }
  if (not cuts.empty() or depth) goto Add;
    // Only if no cuts, no branching
  subgr->Init(lpd,pr->Dim()); // + RHS m.b. set
  do {
    if (fUseLagrange) 
      dbg_outn_(2,"Lagr value: "<<pr->GetLocalLagrValue());
    if (clb < 1e50 && GetLocalLagrBnd() >= clb) {
      if (fUseLagrange && 0==depth)  { // only global yet
//        log_n_(2," lgr="<<GetLocalLagrBnd());
        assertm(lLb <= clb,
          "node:"<<theNode->no<<" localLagr="<<lLv
          <<" > current LPv="<<clv); // not ok yet
        fLPOpt = false; 
        goto RetFalse;
      }
    }
    if (++k > subgr->iterMax) break;
    if (!subgr->UpdateMultipliers // current lgv & lpb! ??
      (pr->GetCurrentLagrValue(),clb,pr->GetBestCol()))
      break; // zero subgradient
  } while (GenCol(subgr->d));
Add:
  if (!AddCols())  // with red. cost < -eps ???
  { dbg_outn_(3,"0C"); goto RetFalse; }
  return true;
RetFalse:
// PricingOver(); // Not here! because no need in PsCosts
  return false;
}
bool BCP::GenCol(d_vec& d) {
  ColSet newCols1;
  double rc;
Again:
  ++ nCG;
  if (depth && clb < gub)
    ++ nCG_TOff;
  if (0 == cntNode % 1000 or 0 == fmod(nCG,1000))
    log_n(3, "\nTailing-off control could save up to "
      << double(nCG_TOff) / nCG << " col. gens.");
  depth ? tmCG2.start() : tmCG1.start();
  try {
    pr->GenCol(newCols1,d);
  } catch (int ) {
    ++ nTooLongColGen;
    DeleteSomeCuts(1);
    SolvePrimal(); // get simpl mult for new cut set
    llvLast = llbLast = -1e100; //to prevent colgen abort
    goto Again;
  }
  depth ? tmCG2.stop() : tmCG1.stop();
  if (0==newCols1.cs.size()) {
//    if (status==none) status = ok; // not the global if.
    return false;
  }
// Now taking only columns with actual negative red.cost
// (optionally) ____ !!!!
  if (not fOriginalRedCost) { // taking all
    for_each_in(newCols1.cs,itc,ColSet::iterator) {
      newCols.cs.insert(*itc);
    }
    // may cause duplications in Soplex.
    // But if good, organize a search tree for ALL cols.
  } else {
    int nGood=0; // But this is false with cuts!!!
    for_each_in(newCols1.cs,itc,ColSet::iterator) {
      ClearCoefs(); // how in BBCuts? No branches cols thr
      rc = pr->CalcRedCost((Column*)&*itc,&cuts,d);
      dbg_outn(5,"Recalc red cost: "<<rc);
      if (rc < - pr->GetRCEps()) {// d !!
        { newCols.cs.insert(*itc); ++ nGood; dbg_outn_(3,"++"); }
      }
    }
    if (0==nGood)
      return false;
  }
  return true;
}
bool BCP::AddCols() {
  bool ret = true;
  int n=lp->NCols();
  if (! AddCols(newCols)) ret = false;
  // See GenCol for checking
  n = lp->NCols() - n;
  dbgcout(5.5,subgr->iterMax<<':'<<newCols.cs.size()<<'/'<<n<<' ');

  if (OUTP_LEV__ > 5) {
    PRINT_LOG("Total N Cols: "<<lp->NCols()<<". Just tried to add these cols:");
    for_each_in(newCols.cs,ic,ColSet::iterator) {
      pr->PrintColumn(mylog,(Column*)&*ic);
      PRINT_LOG__('\n');
    }
  }
  newCols.cs.clear();
  return ret;
}

/// Called when we finish CG, either LP opt or tail off
void BCP::PricingOver() {
  // int i;  // just updates local/global lpv basically
  if (LP::opt == lp->status) {
    UpdateLocalLPV(); // from the current
    UpdateGlobalLPV();
  }
  if (llb != llbLast) iterLLBLast = iter;
    // for tailing off
  llbLast = llb; llvLast = llv;
  if (0==depth && 0==iter) {
    timeLP0 = timer.userTime();
    if (LP::opt == lp->status) {
      status = LPBnd;
      glb0 = glb; glv0 = glv; fLPOpt0 = fLPOpt;
    }
    lpCols = lp->NCols();
    if (fPrintContSol0) {
      GetEnv().GetLog2() 
        << "\n\n\nThe initial LP solution: time=" << timeLP0 << "\n\n";
      PrintContSol();
    }
  }
  dbg_outn_(2,'='<<setprecision(16)<<llv<<setprecision(6));
  dbg_outn_(3, " (Local LP value)");
/*  for (i=0;i<lpx.size();++i) {
    if (not cols[i].slackCut)     // not each get col
    if (lpx[i] > pr->GetXEps() and GetCol(i)->fInfeasible())
      dbg_outn(1, "Infeas column has value "<<lpx[i]);
  }*/
  CheckGivenSolution();
}

void BCP::InitBasis() {
  ColumnList bas;
  pr->HeuristicLPBasis(bas);

  if (OUTP_LEV__ > 4) {
    PRINT_LOG("Basic cols.");
    for_each_in(bas.cl,ic,ColList::iterator) {
      pr->PrintColumn(mylog,&*ic);
      PRINT_LOG__('\n');
    }
    PRINT_LOG__("Constr signs: ");
    for (int i=0;i<pr->validsign.size(); ++i)
      PRINT_LOG__(int(pr->validsign[i]));
    PRINT_LOG__('\n'<<flush);
  }
  lp->InitConstr(pr->b); // and objective
  AddCols(bas);
  AddSlacks();
//  SolvePrimal(); // not nec
}

void BCP::AddSlacks()
{
  int i;
//  lp->InitConstr(pr->b);
  for (i=0;i<pr->Dim(); ++i)
    if (pr->validsign[i]) { // if not an equality
    Column cl;
    cl.MakeSlack(i,
      (pr->validsign[i] > 0) ? -1 : 1 );
    // NOT ADDING TO INDEX ???  DOCH!!
    Column * pc = allcols.Add(cl);
    if (pc) { // there may be valid patterns as slacks
      AddColToLP(pc,colpool.size()); // --the index
      colpool.push_back(pc); // --in the main pool
//      cols.push_back(
  //      ColId(colpool.size()-1, NULL)); // the local
    }
  }
  // ADD INFEASIBLE COLUMNS, NEEDED FOR BRANCHING
  for (i=0;i<pr->Dim(); ++i)
    if (pr->validsign[i] >= 0) { // =b or >=b
    fInfeasColsIn = 1; // !!!!!!! here only
    Column cl;
    cl.MakeSlack(i, 1);
    cl.fInfeas = 1;
    cl.SetObj(maxObjCoef * M__);
    // NOT ADDING TO INDEX ???  DOCH!!
    Column * pc = allcols.Add(cl);
    if (pc) { // there may be valid patterns as slacks
      AddColToLP(pc,colpool.size()); // --the index
      colpool.push_back(pc); // --in the main pool
//      lp->ChangeVarBnds(cols.size()-1,0,0);
  //    log_n_(3,"Fixing infeas. slack for constraint "<<i);
//      cols.push_back(
  //      ColId(colpool.size()-1, NULL)); // the local
    }
  }
}

bool BCP::AddColToLP(Column * c,int j) {
//  c->Sort(); // !!! before

/*if (OUTP_LEV__ >= 5) {
  cout << "addcol: ";
  for (int i=0;i<matind.size();++i)
    cout << matind[i] << ':' << matval[i] << ' ';
  cout << endl;
}*/

  double obj = /*fInfeasColsIn ? 0 // Phase I
    :*/ c->GetObj();
  double lb = 0, // !!! c->GetLB(),
    ub = (fLocalUB && not c->IsSlack()
      && pr->BranchingPossible()) // CSP2
      ? LocalUpperBound(c)  : 1e100;
    // for fLocalUB, could ask Problem additionally
//  int cmatbeg = 0;
//  pr->NotifyAddCol(c); // !!!!! no it knows the set

  if (OUTP_LEV__ >= 2.5) {
    log__(" ADD COL ["<<j<<"]: ")
    for_each_in(c->id,iid,Column::iterator)
      log__(' '<<iid->i<<':'<<iid->d);
    log__(" ub="<<ub);
    log_ln(" obj="<<obj);
  }

  matind.clear(); matval.clear();

  for_each_in(c->id,iid,Column::const_iterator) { // GNU
    matind.push_back(iid->i);
    matval.push_back(iid->d);
  }
//  AddCutCoefs(c,cmatind,cmatval);
  for_each_in(cuts,ic,CutList::iterator) {
    double v=(*ic)->Calc__(c,j);
    if (v) {
      matind.push_back(ic-cuts.begin()+pr->Dim());
      matval.push_back(v);
    }
  }

  int nnz = matind.size();

  lp->AddCol(obj,nnz,
      &matind.front(),&matval.front(),lb,ub);
  cols.push_back(
    ColId(j, NULL)); // the local
// storing column already done.
  return true;
}

void BCP::DelCol(int i) { // IN THE LP NUMERATION
  lp->DelCol(i);
  cols.erase(cols.begin()+i);
}

// Only non-cut slack cols
/// Poss that adding a col with index > rawCoef.size() ?
/// betw nodes with diff cut sets? No 'cause when
/// adding a cut, all coefs are calc-d
bool BCP::AddCols(ColSet &cs) {
  bool fAllExist = true;
  Column *pCol;
  int index; // in the main pool
  assert(colpool.size() == allcols.size());
  for_each_in(cs.cs,ic,ColSet::iterator) {
//    ic->Sort(); must be sorted before adding

    index = colpool.size();
  if (nReducePool) {
    Column * c2 = allcols.Find(*ic);
    if (c2 != NULL) {
//      if (not c2->Hidden())
  //      for_each_in(c2->id,iid,Column::iterator)
    //      log__(' '<<iid->i<<':'<<iid->d);
      //log__('\n');
      assert(c2->Hidden() || c2->IsSlack());
      pCol=c2;
      index = c2->GetIndex();
    }
    else {
      pCol = allcols.Add(*ic);
      assert(pCol);
      pCol->index = index;
    }
    if (pCol) pCol->nHidden = 0; // !!!
  } else
    pCol = allcols.Add(*ic);

    if (pCol) {
      if (colpool.size() == index)
        colpool.push_back(pCol); // --in the main pool
      AddColToLP(pCol,index); // --the index
      if (not pCol->fInfeasible()
        and not pCol->IsSlack())
        if (maxObjCoef < pCol->GetObj()) {
          maxObjCoef = pCol->GetObj();
	}
//      cols.push_back(
  //      ColId(colpool.size()-1, NULL)); // the local
      fAllExist = false;
  // CALC COEFS OF ALL CUTS FOR THIS COL (needed for viol)
      for_each_in(cutpool,ipc,CutRefSet::iterator)
        (*ipc)->Calc__(pCol,index);
    }
  }
  if (OUTP_LEV__ > 5.99) lp->WriteModel("AddCols.lp");
  return !fAllExist;
}

bool BCP::AddCols(ColList &cl) {
  ColSet cs;
  for_each_in(cl.cl,ic,ColList::const_iterator)
    cs.cs.insert(*ic);
  return AddCols(cs);
};

SS_END_NAMESPACE__
