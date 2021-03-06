/* This program implements the algorithm described in Mike Moffitt's
 IJCAI paper. This version runs KK to get the initial solution. In
 this version, SEARCH doesn't return the best solution cost, but
 simple sets the global variable BESTSOFAR to the cost of the best
 solution found so far. This version doesn't pass MAXSOFAR as an
 argument to the search routine, but maintains a global array of its
 value, indexed by k, the number of remaining subsets. */

#include "Moffitt.hpp"

#include <stdio.h>                                    /* standard I/O library */
#include <math.h>                                     /* mathematical library */
#include <iostream>

namespace partition {

uint64_t bestsofar; /* largest subset sum in best solution found so far */
uint64_t maxsofar[MAXK + 1]; /*largest sum of subsets K-index in current solution*/

/* SORT takes an array of longs, and the length of the array, and
 sorts the array in decreasing order using insertion sort. */

void sort(uint64_t a[MAXN], int n)

{
	int i, j; /* indices into array for sorting */
	uint64_t temp; /* temporary value for swapping */

	for (i = 1; i < n; i++) {
		temp = a[i];
		for (j = i - 1; j >= 0; j--)
			if (a[j] < temp)
				a[j + 1] = a[j];
			else
				break;
		a[j + 1] = temp;
	}
}

/* INSERT takes an array index FIRST, a new subset sum VECTOR, and an array A of
 vectors sorted in decreasing order of largest elements, and modifies array A
 by inserting the new vector in sorted order by largest element. */

void insert(int first, int n, int k, uint64_t vector[MAXK], uint64_t vects[MAXN][MAXK])

{
	int i; /* index into array of subpartitions */
	int j; /* index to elements of vector */

	for (i = first; i < n - 1; i++) /* insert new partition in sorted list */
		if (vector[0] < vects[i + 1][0]) /* haven't found correct place yet */
			for (j = 0; j < k; j++) /* for each element of vector */
				vects[i][j] = vects[i + 1][j]; /* copy current partition up in order */
		else
			break; /* found correct place, exit loop */

	for (j = 0; j < k; j++) /* for each element of vector */
		vects[i][j] = vector[j];
} /* insert new subpartition in order */

/* GREEDY takes an array A of integers, sorted in decreasing order,
 the number N of integers in the array, and the number of subsets K,
 and returns the largest subset sum in a greedy solution to the
 K-way partitioning of the integers. */

uint64_t greedy(uint64_t a[MAXN], int n, int k)

{
	int i, j; /* utility indices */
	uint64_t newsum; /* sum of adding next number to smallest subset */
	uint64_t set[MAXK]; /* subsets constructed */

	for (i = 0; i < k; i++) /* put first k integers in subsets */
		set[i] = a[i];
	for (i = k; i < n; i++) /* for each remaining integer */
	{
		newsum = set[k - 1] + a[i]; /* put number in smallest subset */
		for (j = k - 2; j >= 0; j--)
			if (set[j] < newsum)
				set[j + 1] = set[j];
			else
				break;
		set[j + 1] = newsum;
	}
	return (set[0]);
}

/* SEARCH is the main recursive routine.  It is used both to iterate
 through successive subsets and to perform the inclusion-exclusion
 search for a given subset.  It takes as arguments the number of
 sets K in which to divide the remaining numbers, the sum of the
 numbers included so far in the current subset (SUBSUM), an array of
 the numbers excluded from the current subset (EX), the number of
 excluded numbers (NUMEX), the sum of the excluded numbers (SUMEX),
 an array A of remaining integers to be included or excluded, sorted
 in decreasing order, the length N of the array, the SUM of all the
 numbers in the array, and a LOWER bound on the subset sum currently
 being constructed.  The upper bound on any subset sum is the global
 variable BESTSOFAR.  As a side effect, it sets BESTSOFAR to the
 largest subset sum in the best complete solution found so far. */

void search(int k, /* number of subsets in which to divide remaining numbers */
uint64_t subsum, /*sum of integers included so far in current subset */
uint64_t ex[MAXN], /* array of integers excluded from this subset */
int numex, /* number of excluded numbers */
uint64_t sumex, /* sum of excluded numbers */
uint64_t a[MAXN], /* integers not assigned to any subsets yet */
int n, /* length of integer array A */
uint64_t sum, /* sum of integers in A */
uint64_t lower) /* the minimum sum required for current subset */

{
//	cout << "K: " << k << " ss: " << subsum << " a[0]: " << a[0] << " bsf: " << bestsofar << endl;
	uint64_t newsum; /* maximum subset sum in just completed solution */
	uint64_t newlower; /* new lower bound for next subset sum */
	uint64_t newex[MAXN]; /* new array for excluded numbers */
	uint64_t last; /* sum of last subset sum when greedy is optimal */
	if (n == k && numex == 0 && k > 2) /* greedy is optimal at this point */
	{
		last = a[n - 2] + a[n - 1]; /* smallest two numbers go in last subset */
		if (last >= subsum)
			newsum = last; /* last subset is larger than first */
		else
			newsum = subsum; /* first subset is larger than last */
		if (newsum <= maxsofar[k])
			bestsofar = maxsofar[k]; /* previous subset sums larger */
		else if (newsum < bestsofar)
			bestsofar = newsum; /* new sol better than best so far*/
		return;
	} /* end search */

	if (n == 1) /* only one number left to be considered */
	{
		if (k == 2) /* completing penultimate subset */
		{
//			cout << "k==2: ss = " << subsum << " ss_ex = " << sumex << " a[0]: " << a[0] << endl;
			if (subsum <= sumex) {/* current subset is smaller than complement */
				if (subsum + a[0] >= sumex)
					newsum = subsum + a[0]; /* place in smaller set */
				else
					newsum = sumex; /* complement has larger sum */
			} else {/* add last number to complement set */
				if (sumex + a[0] >= subsum)
					newsum = sumex + a[0]; /* new set is larger */
				else
					newsum = subsum; /* current subst is larger */
			}

			if (newsum <= maxsofar[k]) {
				bestsofar = maxsofar[k]; /*previous subset sums larger*/
			} else if (newsum < bestsofar) {
				bestsofar = newsum; /*better solution found*/
			}
//			cout << "bsf: " << bestsofar << endl << "---------------------------" << endl;
			return;
		} /* return best solution found */

		/* one number left to assign, but not penultimate subset */
		newsum = subsum + a[0]; /* include last number in current subset */
		if (newsum < bestsofar && newsum + ex[numex - 1] > maxsofar[k]) /* not dominated */
		{
			if (newsum > maxsofar[k])
				maxsofar[k - 1] = newsum; /*maxsum of completed subsets*/
			else
				maxsofar[k - 1] = maxsofar[k];
			newlower = sumex - (k - 2) * (bestsofar - 1); /*lb on next subset */
			if (newlower < bestsofar) /*lower bound must be better than best solution so far*/
			{
//				// *** DEBUGGING BEGIN ***
//				cout << "Recurse " << k-1 << endl;
//				cout << "LB: " << lower << " BSF: " << bestsofar << endl;
//				for (int i = 0; i < n; i++) {
//					cout << ex[i] << " ";
//				}
//				cout << "sum: " << subsum << endl;
//				cout << "EX1: ";
//				for (int i = 0; i <= numex; i++) {
//					cout << ex[i] << " ";
//				}
//				cout << endl << endl;
//				// *** DEBUGGING END ***
				search(k - 1, ex[0], newex, 0, 0ll, &ex[1], numex - 1, sumex - ex[0], newlower);
			}
			if (bestsofar <= maxsofar[k])
				return;
		} /* no need to continue searching*/

		if (newsum >= bestsofar && subsum >= lower) /*lb reached, inclusion doesn't dominate*/
		{
			ex[numex] = a[0]; /* move last number to excluded list */
			if (subsum > maxsofar[k])
				maxsofar[k - 1] = subsum; /*new maximum sum of completed subsets*/
			else
				maxsofar[k - 1] = maxsofar[k];
			newlower = sumex + a[0] - (k - 2) * (bestsofar - 1); /*lb on next subset*/
			if (newlower < bestsofar) {
				// // *** DEBUGGING BEGIN ***
				// cout << "Recurse " << k-1 << endl;
				// cout << "LB: " << lower << " BSF: " << bestsofar << endl;
				// for (int i = 0; i < n; i++) {
				// 	cout << ex[i] << " ";
				// }
//				 cout << "sum: " << subsum << endl;
//				 cout << "EX2: ";
//				 for (int i = 0; i <= numex; i++) {
//				 	cout << ex[i] << " ";
//				 }
//				 cout << endl << endl;
				// // *** DEBUGGING END ***
				search(k - 1, ex[0], newex, 0, 0ll, &ex[1], numex, sumex + a[0] - ex[0], newlower);
			}
		}
		return;
	}
	/* more than one number left to assign */

	if (subsum + a[0] < bestsofar) /* including next number doesn't exceed upper bound */
	{
		search(k, subsum + a[0], ex, numex, sumex, &a[1], n - 1, sum - a[0], lower);
		if (bestsofar <= maxsofar[k])
			return; /* solution <= max of previous subsets*/

		// if (subsum + a[0] < bestsofar && subsum + a[0] >= lower) /* Moffitts pruning rule */
		// 	lower = subsum + a[0] + 1;
	}

	if (subsum + sum - a[0] >= lower) /* if next number is excluded, can still make lb */
	{
		ex[numex] = a[0]; /* add next number to excluded array */
		search(k, subsum, ex, numex + 1, sumex + a[0], &a[1], n - 1, sum - a[0], lower);
	}
}

uint64_t executeMoffitt(const PartitionProblem &problem, ProblemStats &stats) {

	uint64_t S[MAXN]; /* array of integers for each problem instance */
	uint64_t ex[MAXN]; /* array for excluded numbers */
	uint64_t lower; /* lower bound on first subset sum */

	int K = problem.K;                     	// total number of sets to partition numbers into
	int N = problem.N;                    	// problem size: number of values in number array
	uint64_t sumall = problem.sum;         	// sum of all the numbers in the problem instance

	memcpy(S, problem.S, N * sizeof(uint64_t));

	sort(S, N); /* sort numbers in decreasing order */

	bestsofar = kk(S, N, K, sumall); /* compute KK approximation */

	lower = sumall - (K - 1) * (bestsofar - 1); /* lower bound on first subset*/
	//cout << "KK: " << bestsofar << endl;
	//cout << "LB: " << lower << endl;
	//	cout << "Sum: " << sumall << endl;
	maxsofar[K] = 0ll; /* initially, no completed subsets */

	//	void
//	search (int k,      /* number of subsets in which to divide remaining numbers */
//	        uint64_t subsum, /*sum of integers included so far in current subset */
//	        uint64_t ex[MAXN],    /* array of integers excluded from this subset */
//	        int numex,                              /* number of excluded numbers */
//	        uint64_t sumex,                           /* sum of excluded numbers */
//	        uint64_t a[MAXN],        /* integers not assigned to any subsets yet */
//	        int n,                                   /* length of integer array A */
//	        uint64_t sum,                                /* sum of integers in A */
//	        uint64_t lower)       /* the minimum sum required for current subset */
//

	search(K, S[0], ex, 0, 0, &S[1], N - 1, sumall - S[0], lower);

	return bestsofar;
} // end executeMoffitt
} // end namespace
