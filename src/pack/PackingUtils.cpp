/*
 * BinPackingUtils.cpp
 *
 *  Created on: Aug 2, 2012
 *      Author: ethan
 */

#include "PackingUtils.hpp"

using std::cout;
using std::endl;


// *******************
// PackingProblem
// *******************
PackingProblem::PackingProblem(int N0, uint64_t S0[], uint64_t sum0, string problemName0 )
: N(N0), sum(sum0), problemName(problemName0) {
  S = new uint64_t[N];
  memcpy(S,S0,sizeof(uint64_t) * N);
}

PackingProblem::~PackingProblem() {
  delete [] S;
}
PackingProblem::PackingProblem(const string &filename) {
  load(filename);
}


PackingProblem::PackingProblem(const PackingProblem &problem)
: N(problem.N), sum(problem.sum), problemName(problem.problemName) {

  S           = new uint64_t[N];
  memcpy(S,problem.S,N * sizeof(uint64_t));
}

void PackingProblem::readElements(std::ifstream &inFile) {

  S = new uint64_t[N];     // Allocate memory for S
  for (int i=0;i<N;i++) { // Fill in S values
    if (!inFile) {
      cout << "Not enough entries, expected " << N << ", got " << i << endl;
      exit(0);
    }
    inFile >> S[i];
  }

  uint64_t next;     // Make sure no more
  inFile >> next;

  if (inFile) {
    cout << "Entries left over, perhaps you should have used a different PackingProblem struct?" << endl;
    cout << "Next: " << next << endl;
    cout << "Sum : " << sum << endl;
    for (int i=0;i<N;i++) {
      cout << "--> " << S[i] << endl;
    }
    exit(0);
  }

  sum = std::accumulate(S,S + N,(uint64_t) 0);   // sum is sum of all elements in S

}

void PackingProblem::load(string filename) {
  std::ifstream inFile(filename);
  std::getline(inFile,problemName);   // Read problem Name

  inFile >> N >> sum;                       // [N] [sum]
  readElements(inFile);
  inFile.close();
}



// *****************
// BinPackingProblem
// *****************
BinPackingProblem::BinPackingProblem(const string &filename, int lowerBound) {

  //cout << "Reading " << filename << endl;

  std::ifstream inFile(filename);
  std::getline(inFile,problemName);

  inFile >> capacity >> N >> sum;
  readElements(inFile);
}

BinPackingProblem::BinPackingProblem(const PackingProblem &packingProblem, uint64_t capacity, int lowerBound)
: PackingProblem(packingProblem), capacity(capacity), lowerBound(lowerBound)
{}

BinPackingProblem::BinPackingProblem(const uint64_t capacity, int lowerBound) : capacity(capacity), lowerBound(lowerBound) {
  S = NULL;
  N = 0;
  sum = 0;
}

void BinPackingProblem::reset(uint64_t S0[MAXN],size_t N0) {
  if (S != NULL) {
    delete [] S;
  }
  sum = 0;
  N = N0;
  S = new uint64_t[N];

  for (int i=0;i<N;i++) {
    S[i] = S0[i];
    sum += S[i];
  }
}
string BinPackingProblem::toString() const {
  std::ostringstream out;
  out << "Problem Name: " << problemName << endl
      << "Bin Capacity: " << capacity<< endl
      << "# Items     : " << N << endl
      << "Sum         : " << sum << endl;
  return out.str();
}

// *****************************
// Bin Packing Utility Functions
// *****************************

// Upper Bound
uint64_t BFD(uint64_t items[],int solution[],int numItems,uint64_t capacity) {                                /* array of numbers */

  int numBins = 0;  // number of bins used
  int *binum = new int[numItems];  // the original # of the bin in each current position
  uint64_t * binSum = new uint64_t[numItems]; // the sums in the bins so far

  for (int i = 0; i < numItems; i++) {      // for each bin
    binum[i] = i;                           // its original index
    binSum[i] = 0;                          // initialize bins to empty
  }

  for (int i = 0; i < numItems; i++) {          // for each number
    int j = 0;                                  // start with first bin
    while (items[i] + binSum[j] > capacity) {
      j++;                                      // try next if it doesn't fit
    }
    uint64_t newbinsum = binSum[j] + items[i];       // insert into first bin it fits, or last bin
    solution[i] = binum[j];                     // store assignments in solution vector
    if (j >= numBins) {
      numBins++;                                // keep track of how many bins used
    }

    while (j > 0 && binSum[j-1] < newbinsum) {  // resort bins from fullest to emptyest
      binSum[j] = binSum[j-1];                  // element is less full than newbin
      binum[j] = binum[j-1];                    // keep track of original indices
      j--;                                      // see if greater than next bin as well
    }
    binSum[j] = newbinsum;                      // insert new bin in sorted order
    binum[j] = solution[i];                     // keep track of original index of bin
  }

  delete [] binum;
  delete [] binSum;
  return (numBins);
}

uint64_t BFD(uint64_t items[],int numItems,uint64_t capacity) {
  int * solution = new int[numItems];
  uint64_t result = BFD(items,solution,numItems,capacity);
  delete [] solution;

  return result;
}


// Lower bound
uint64_t L2(const uint64_t S[MAXN],         // The elements to be assigned
       int N,               // The number of elements
       uint64_t capacity,        // The capacity of a bin
       uint64_t sum) {           // Sum of all of the elements.

  int first     = 0;        // the index of the largest element
  int last      = N-1;      // the index of the smallest element
  uint64_t waste     = 0;        // the amount of wasted space computed so far
  uint64_t carryover = 0;        // the sum of the elements that haven't been used so far

  while (first <= last) {                     // while there are still elements remaining
    uint64_t empty = capacity - S[first];              // empty space in bin with largest element
    while (last > first && S[last] <= empty){ // elements left that fit
      carryover += S[last];                   // add to total that could fit in this bin
      last--;                                 // try next larger element
    }

    if (carryover < empty) {                  // all that can fit does
      waste += empty - carryover;             // rest is wasted space
      carryover = 0;                          // all carryover is consumed by this bin
    } else {
      carryover -= empty;                     // reduce carryover by as much as fits in bin
    }
    first++;                                  // try next largest element
  }

  uint64_t total = sum + waste;                    // total capacity needed
  uint64_t numBins = total / capacity;             // number of bins needed for this total
  if(numBins * capacity < total) numBins++;   // didn't divide evenly, round up

  return numBins;                             // return lower bound on bins needed
}





