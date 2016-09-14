/*
 * BinPackingTest.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: ethan
 */

#include "BinCompletion.hpp"

// Was playing around with beating Belov for hard28. Possible?

struct Subset {
private :
  std::vector<uint64_t> m_values;
public :
  Subset(const uint64_t values[], size_t N) {
    for (size_t i=0;i<N;i++) {
      m_values.push_back(values[i]);
    }
  }
  void push_back(uint64_t value) {
    m_values.push_back(value);
  }
  string toString() {
    std::ostringstream out;
    out << "{";
    uint64_t sum=0;
    for (size_t i=0;i<m_values.size();i++) {
      if (i != 0) {
        out << ", ";
      }
      out << std::setw(3) << m_values[i];
      sum += m_values[i];
    }
    out << "}  Sum = " << sum ;

    return out.str();
  }

};
void print(uint64_t sum, int numBins, uint64_t capacity) {
  cout << numBins << "      " << (double) sum / (double) numBins << "     " << capacity*numBins - sum<< endl;
}

// x-> Number of bins with 2 elements
// y-> Number of bins with 3 elements
// 2x + 3y = N
// x+y=numBins

// x = numBins-y
// 2(numBins-y) + 3y = N
// y = N - 2*numBins

size_t numXY(int N, int numBins) {
  uint64_t y = N-2*numBins;
  uint64_t x = numBins-y;

  cout << "Num Bins      : " << numBins << endl
       << "2 element Bins: " << x << endl
       << "3 element Bins: " << y << endl << endl;
  return x;
}

void getUndominated2( vector<Subset> &subsets, const uint64_t S[], const size_t N, const uint64_t C) {

  for (size_t i=0;i<N;i++) {
    for (size_t j=i+1;j<N;j++) {
      if (S[i] + S[j] <= C) {
        uint64_t arr[2] = {S[i],S[j]};
        subsets.push_back(Subset(arr,2));
        break;
      }
    }
  }
}


size_t printPairs( const uint64_t S[], const size_t N, const uint64_t C, const size_t total) {

  cout << "----- Capacity: " << C << " -----" << endl;
  int i=0;
  int j=N-1;

  size_t count=0;
  while (i < j) {
    uint64_t sum = S[i] + S[j];
    if (sum == C) {
      cout << S[i] << " + " << S[j] << " = " << C << endl;
      count++;
      i++;
      j--;
    } else if (sum > C) {
      i++;
    } else {
      j--;
    }
  }

  cout << "Count: " << count << endl;
  cout << "Total: " << count+total << endl;
  cout << endl;
  return count;
}


size_t removePerfectPairs(uint64_t *&S, size_t &N, uint64_t C, size_t total) {
  cout << "----- Capacity: " << C << " -----" << endl;

  vector<uint64_t> firstHalf;
  deque<uint64_t> secondHalf;

  int i=0;
  int j=N-1;

  size_t count=0;
  while (i < j) {
    uint64_t sum = S[i] + S[j];
    if (sum == C) {
      cout << S[i] << " + " << S[j] << " = " << C << endl;
      count++;
      i++;
      j--;
    } else if (sum > C) {
      firstHalf.push_back(S[i]);
      i++;
    } else {
      secondHalf.push_front(S[j]);
      j--;
    }
  }
  if (i == j) {
    firstHalf.push_back(S[i]);
  }
  firstHalf.insert(firstHalf.end(),secondHalf.begin(),secondHalf.end());

  N = firstHalf.size();
  memcpy(S,&firstHalf[0],N * sizeof(uint64_t));


  cout << "Count: " << count << endl;
  cout << "Total: " << count+total << endl;
  cout << endl;
  return count;
}

void executeTest( const BinPackingProblem &problem,ProblemStats &stats, const PackingOptions &options) {

  uint64_t C = problem.capacity;
  size_t N = problem.N;
  uint64_t *S = new uint64_t[N];
  uint64_t lb = ( problem.sum  + problem.capacity-1) /  problem.capacity;
  uint64_t maxWaste = C *lb - problem.sum;


  memcpy(S,problem.S,N * sizeof(uint64_t));


  vector<Subset> subsets;
  getUndominated2(subsets,S,N,C);
  cout << endl;

  cout << "# Bins  Value/Bin   Waste" << endl;
  print(problem.sum,lb,C);

  cout << endl;
  //size_t num2 = numXY(problem.N,lb);


  size_t total=0;
//  total += printPairs(S,N,C, total);
//  total += printPairs(S,N,C-1, total);
//  total += printPairs(S,N,C-2, total);
//  total += printPairs(S,N,C-3, total);
  int wastePerPair=0;
  int waste=0;
  while (waste < (int) maxWaste) {

    size_t count = removePerfectPairs(S,N,C,total);
    total += count;
    C--;
    waste += (wastePerPair * count);
    wastePerPair++;

  }

  cout << "Waste: " << waste << endl;

  for (size_t i=0;i<subsets.size();i++) {
    cout << subsets[i].toString() << endl;
  }

}
