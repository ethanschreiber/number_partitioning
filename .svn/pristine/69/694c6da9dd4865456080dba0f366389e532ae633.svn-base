/* clear && g++ generateAll.cpp && ./a.out 4 2 5 9 23 47 51 54 63 > a.txt */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <algorithm>

#include <math.h>
using namespace std;

typedef vector<int> IntVec;
typedef vector<IntVec> IntVecVec;

const int SUM_WIDTH = 3;
const int ELEMENT_WIDTH = 2;

struct Comparator {
	
	bool operator()(const IntVec &S1, const IntVec &S2) {
		int sum1 = std::accumulate(S1.begin(),S1.end(),0);
		int sum2 = std::accumulate(S2.begin(),S2.end(),0);
		return sum1 < sum2;
	}
};

void generateAll(IntVec &S, IntVecVec &result) {
  int numSets = 1 << S.size();
	cout << "Num Sets: " << numSets << endl;

	for (int i=0;i<numSets;i++) {
		IntVec set;
		for (int j=0;j<S.size();j++) {
			if ((i >> j) & 1) {
				set.push_back(S[j]);
			}
		}
		result.push_back(set);
	}
	std::sort(result.begin(),result.end(),Comparator());
}

int main(int argc, char *argv[]) {

	if (argc < 3) {
		cout << "Usage: " << argv[0] << "k S1 S2 S3 ..." << endl;
		exit(0);
	}
	IntVec S;

	int K = atoi(argv[1]);


	for (int i=2;i<argc;i++) {
		S.push_back(atoi(argv[i]));
	}
	const int totalSum = (int) std::accumulate(S.begin(),S.end(), 0);
	const int perfect = (int) ceil((double) totalSum / (double) K);
	const int ALL_SETS_WIDTH = S.size() * (ELEMENT_WIDTH+1);

  cout << "K       : " << K << endl;
  cout << "Sum     : " << totalSum << endl;
  cout << "Perfect : " << perfect << endl;
	cout << "Input   : {";
	for (size_t i=0;i<S.size();i++) {
		if (i != 0) { cout << ", "; }
		cout << S[i];
	} cout << "}" << endl;


	IntVecVec result;
	generateAll(S,result);

	cout << endl
	     << " SUM | Set" << endl
	      << " " << string(SUM_WIDTH+1,'-') << "+" << string(ALL_SETS_WIDTH+1,'-') << endl;
	int idx = 0;
	int perfectIdx = -1;
	for (size_t i=0; i <result.size();i++) {
		IntVec &set = result[i];
		int sum = std::accumulate(set.begin(),set.end(),0);
		cout << " " << std::right << std::setw(SUM_WIDTH) << sum << " | ";

		std::ostringstream out;
		for (size_t j=0;j<set.size();j++) {
			if (j != 0) { out << " "; }
			out << std::setw(ELEMENT_WIDTH) << set[j];
		}
		cout << std::left << std::setw(ALL_SETS_WIDTH) << out.str() << "| ";

		if (sum > perfect) {
			if (perfectIdx == -1) {
				perfectIdx = i;
			}
		  idx++;
		  int CMin = totalSum - (K-1) * (sum);
		  cout << std::setw(2) << idx << " (CMin = " << CMin <<")";
		  if (CMin < 0) {
		  	break;
		  }

		}
		cout << endl;
	}


	cout << endl << endl;
	cout << "Iterative Weakining" << endl
			 << "-------------------" << endl << endl;
	size_t minIdx=perfectIdx-1;
	IntVec *minPtr = &result[minIdx];
	int minSum = std::accumulate(minPtr->begin(),minPtr->end(),0);
	for (size_t i=perfectIdx; i <result.size();i++) {
		IntVec &set = result[i];
		int sum = std::accumulate(set.begin(),set.end(),0);
		cout << " " << std::right << std::setw(SUM_WIDTH) << sum << " | ";

		std::ostringstream out;
		for (size_t j=0;j<set.size();j++) {
			if (j != 0) { out << " "; }
			out << std::setw(ELEMENT_WIDTH) << set[j];
		}
		cout << std::left << std::setw(ALL_SETS_WIDTH) << out.str() << "| ";


		if (perfectIdx == -1) {
			perfectIdx = i;
		}
		idx++;
		int CMin = totalSum - (K-1) * (sum);

		if (CMin < 0) {
			break;
		}


		cout << " (CMin = " << std::setw(2) << CMin <<"): ";
		while (minIdx >= 0 && minSum > CMin) {
			std::ostringstream out;
			for (size_t j=0;j<result[minIdx].size();j++) {
				if (j != 0) { out << " "; }
				out << std::setw(ELEMENT_WIDTH) << result[minIdx][j];
			}
			cout << "{" << out.str() << " } ";

			minIdx--;
			minPtr = &result[minIdx];
			minSum = std::accumulate(minPtr->begin(),minPtr->end(),0);
		}



		cout << endl;
	}
}
