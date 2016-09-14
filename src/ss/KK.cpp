/*
 * kk.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: ethan
 */

#include "KK.hpp"
#include <iomanip>
using std::cout;
using std::endl;


// ============================================================================
// INSERT takes an array index FIRST, a new subset sum VECTOR, and an array A of
// vectors sorted in decreasing order of largest elements, and modifies array A
// by inserting the new vector in sorted order by largest element.
// ============================================================================
void kkInsert(int first, int n, int k, uint64_t vector[],	uint64_t **vects) {
	int i;                                  // index into array of subpartitions

	for (i = first; i < n - 1; i++) {       // insert new partition in sorted list
		if (vector[0] < vects[i + 1][0]) {    // haven't found correct place yet
			for (int j = 0; j < k; j++) {       // for each element of vector
				vects[i][j] = vects[i + 1][j];    // copy current partition up in order
			}
		} else {
			break;                              // found correct place, exit loop
		}
	}
	for (int j = 0; j < k; j++) {           // for each element of vector
		vects[i][j] = vector[j];              // insert new subpartition in order
	}
}

// ============================================================================
// KK takes an array A of integers, sorted in decreasing order, its length N,
// the number of partitions K, and the sum of all the numbers, and returns the
// maximum subset sum in the KK approximation of the best K-way partition.
// ============================================================================
uint64_t kk(const uint64_t S[], int n, int k, uint64_t sum) {
	uint64_t **vectorArray;           							// array of vectors N long and K wide
	vectorArray = new uint64_t*[n];									// Allocate rows of
	for (int i=0;i<n;i++) {
		vectorArray[i] = new uint64_t[k];
	}
	uint64_t *vectorCombo = new uint64_t[k]; 			// combination of first two vectors


	if (k >= n) {
		return (S[0]);             							// largest subset sum is largest number
	}

	for (int i = 0; i < n; i++) {             // create initial array of vectors
		vectorArray[i][0] = S[i];                			// copy to vectors in decreasing order
		for (int j = 1; j < k; j++) {
			vectorArray[i][j] = 0;
		}
	}

	for (int i = 0; i < n - 1; i++) {

//    for (int idx=i;idx<n;idx++) {
//      cout << std::setw(2) << idx << ": ";
//      for (int j=0;j<k;j++) {
//        cout << std::setw(2) << vectorArray[idx][j] << " ";
//      }
//      cout << endl;
//    }cout << endl;

		for (int j = 0; j < k; j++) {           // for each subset sum
			vectorCombo[j] = vectorArray[i][j] +
									vectorArray[i + 1][k - j - 1];	// combine particular subsets
		}
		std::sort(vectorCombo, vectorCombo + k,
							std::greater<uint64_t>()); 		// Sort in ascending order
		//sort (vector, k);                   	// sort the vector in decreasing order
		for (int j = 0; j < k; j++) {             	// for each subset sum
			vectorCombo[j] = vectorCombo[j] -
									vectorCombo[k - 1];     				// subtract smallest subset sum
		}
		kkInsert(i + 1, n, k, vectorCombo, vectorArray);   	// insert new vector in sorted order
	}

//	cout << "Final Vec: ";
//	for (int i=0;i<k;i++) {
//	  cout << std::setw(2) << vectorCombo[i] << " ";
//	} cout << endl;

	for (int j = 0; j < k; j++) {       			// compute largest subset sum from differences
		sum -= vectorCombo[j];                	// subtract relative subset sums
	}

	uint64_t max = (sum / k) + vectorCombo[0];            	// absolute maximum subset sum

//	cout << "Sum: " << sum << endl;
//	cout << "CKK: " << max << endl;
	// Free memory
	for (int i=0;i<n;i++) {
		delete [] vectorArray[i];
	}
	delete [] vectorArray;
	delete [] vectorCombo;
	return max;
}

//// ============================================================================
//// M a i n
//// ============================================================================
//int main() {
//	const int N = 50;
//	const int K = 10;
//	uint64_t a[N] = { 269953206533712, 264767189672722, 260109545946063,
//			246810908606432, 245982129682025, 245014179504882, 244241212467227,
//			237856245670097, 222213478745249, 221182828081148, 218646933430302,
//			217029100824126, 216761539341028, 211078642492280, 209725625483456,
//			206185946358013, 201406209410675, 198494903516273, 194835333224023,
//			162496491130133, 162069970321160, 160042568757305, 153468045305767,
//			148994111536649, 128303367473637, 125550761873867, 122701314126126,
//			121206344577643, 110983831113855, 109968493374502, 107574789254391,
//			103798477237526, 99565637727610, 88269004309149, 77548868773837,
//			74282998913434, 73285133115476, 73245981694404, 61541592842792,
//			48083817484545, 45933763597706, 43554508934966, 43406518304028,
//			35920996834968, 30947276960908, 30227846306591, 27126209522211,
//			23363618230530, 23138843563715, 7198701771445 };
//
//	uint64_t sum = std::accumulate(a, a + N, (uint64_t) 0);
//
//	cout << "KK: " << kk(a, N, K, sum) << endl;
//}
