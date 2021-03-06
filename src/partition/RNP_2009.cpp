





/* This program optimally solves arbitrary K-way partitioning problems while
   minimizing the maximum subset sum objective function.  It is intended to
   generalize rnp4max.c and rnp5max.c from the IJCAI09 paper to arbitrary k-way
   partitioning.  First, the KK heuristic is used to compute an approximate
   K-way partition.  Then, if K is even, all the numbers are partitioned
   approximately in half, and then each half is recursively partitioned K/2
   ways.  Alternatively, if K is odd, a first subset if computed, and then the
   remaining elements are recursively partitioned K-1 ways.  The bounds on an
   individual subset sum range from an upper bound equal to sum/k, down to a
   lower bound large enough that the remaining elements can be partitioned into
   a set of subsets each of whose sum is less than minmax-1. */

#include <stdio.h>                                    /* standard I/O library */
#include <math.h>                                     /* mathematical library */
#include <stdio.h>
#include "RNP.hpp"
namespace partition {

#define MAXK 12         /* maximum number of subsets to partition values into */
#define MAX_RNP2009 60                 /* maximum number of values being partitioned */


long random(void);                                 /* random number generator */
void srandom(unsigned int seed);   /* function to see random number generator */

int K;                            /* number of sets to partition numbers into */
int N;                      /* problem size: number of values in number array */
int TRIALS;                                        /* number of trials to run */

long long sorted[MAX_RNP2009];        /* original numbers sorted in decreasing order */
long long sumall;                              /* sum of all original numbers */
long long minmax;         /* largest subset sum in best solution found so far */
long long nodes[MAXK];                     /* number of K-way partition calls */
long long A[MAX_RNP2009];        /* global array of numbers to be partitioned by CKK */
long long alpha;      /* global variable for smallest difference found by CKK */


/* SORTDOWN takes an array of integers, and the length of the array, and sorts
   the array in decreasing order using selection sort. */

void sortdown (long long a[MAX_RNP2009], int n)

{int i,j;                                   /* indices into array for sorting */
 long long temp;                              /* temporary value for swapping */

 for (i = 0; i < n-1; i++)                        /* pointer to first element */
   for (j = i+1; j < n; j++)                     /* pointer to second element */
     if (a[i] < a[j])                            /* elements are out of order */
       {temp = a[i];
        a[i] = a[j];                                             /* swap them */
        a[j] = temp;}}

/* INSERT takes an array index FIRST, a new subset sum VECTOR, and an array A of
   vectors sorted in decreasing order of largest elements, and modifies array A
   by inserting the new vector in sorted order by largest element. */

void insert (int first, int n, int k, long long vector[MAXK], long long vects[MAX_RNP2009][MAXK])

{int i;                                  /* index into array of subpartitions */
 int j;                                        /* index to elements of vector */

 for (i = first; i < n-1; i++)         /* insert new partition in sorted list */
   if (vector[0] < vects[i+1][0])          /* haven't found correct place yet */
     for (j = 0; j < k; j++)                    /* for each element of vector */
       vects[i][j] = vects[i+1][j];     /* copy current partition up in order */
   else break;                              /* found correct place, exit loop */

 for (j = 0; j < k; j++)                        /* for each element of vector */
   vects[i][j] = vector[j];}              /* insert new subpartition in order */

/* KK takes an array A of integers, sorted in decreasing order, its length N,
   the number of partitions K, and the sum of all the numbers, and returns the
   maximum subset sum in the KK approximation of the best K-way partition. */

long long kk (long long a[MAX_RNP2009], int n, int k, long long sum)

{long long vects[MAX_RNP2009][MAXK];           /* array of vectors N long and K wide */
 long long vector[MAXK];                  /* combination of first two vectors */
 int i, j;                                                 /* utility indices */
 long long max;                       /* absolute value of maximum subset sum */

 if (k >= n) return (a[0]);           /* largest subset sum is largest number */

 for (i = 0; i < n; i++)                   /* create initial array of vectors */
   {vects[i][0] = a[i];                  /* copy to vects in decreasing order */
   for (j = 1; j < k; j++)
     vects[i][j] = 0;}

 for (i = 0; i < n-1; i++)
   {for (j = 0; j < k; j++)                            /* for each subset sum */
     vector[j] = vects[i][j] + vects[i+1][k-j-1];/* combine particular subsets*/
   sortdown (vector, k);               /* sort the vector in decreasing order */
   for (j = 0; j < k; j++)                             /* for each subset sum */
     vector[j] = vector[j] - vector[k-1];     /* subtract smallest subset sum */
   insert (i+1, n, k, vector, vects);}   /* insert new vector in sorted order */

 for (j = 0; j < k; j++)       /* compute largest subset sum from differences */
   sum = sum - vector[j];                    /* subtract relative subset sums */
 max = (sum / k) + vector[0];                  /* absolute maximum subset sum */
 return (max);}

/* CKK takes an array of long long numbers A, its length N, and their TOTAL sum,
   and finds the best partition, leaving the resulting difference in the global
   variable ALPHA. It runs branch-and-bound, starting with the Karmarkar-Karp
   solution. At each point, the largest two remaining numbers are selected, and
   replaced with either their difference or their sum, representing assigning
   them to different sets or the same set, respectively. */

void ckk (
long long a[MAX_RNP2009],                                        /* array of numbers */
int n,                                         /* number of elements in array */
long long total)                                        /* sum of all numbers */

{long long diff2;                          /* difference of largest 2 numbers */
 long long sum2;                                  /* sum of largest 2 numbers */
 long long difference;                             /* difference of partition */
 long long b[MAX_RNP2009];                   /* new copy of list for recursive calls */
 long long rest;                        /* sum of all elements except first 2 */
 int i;                                                   /* index into array */

 diff2 = a[0] - a[1];                   /* difference of two largest elements */
 for (i = 2; i < n; i++)          /* copy list and insert difference in order */
   if (diff2 < a[i])
     b[i-2] = a[i];                       /* new number is less, keep copying */
   else break;                              /* found correct place, exit loop */
 b[i-2] = diff2;                             /* insert new element into array */
 for (i = i; i < n; i++)                           /* copy remaining elements */
   b[i-1] = a[i];

 if (n == 5)                   /* when only 4 numbers are left, KK is optimal */
   {diff2 = b[0] - b[1];                  /* difference of 2 largest elements */
    if (diff2 < b[2]) difference = b[2] - b[3] - diff2;    /* b[2] is biggest */
    else difference = diff2 - b[2] - b[3];                /* diff2 is biggest */
    if (difference < 0) difference = -difference;      /* take absolute value */
    if (difference < alpha)            /* better than best previous partition */
      alpha = difference;}                     /* reset alpha to better value */

 else
   {rest = total - a[1] - a[1] - b[0];    /* sum of all elements except first */
    if (b[0] >= rest)                    /* if largest element >= sum of rest */
      {if (b[0] - rest < alpha)               /* best better than best so far */
        alpha = b[0] - rest;                   /* reset alpha to better value */
      return;}            /* if difference is larger than rest, so is the sum */
    ckk (b, n-1, total - a[1] - a[1]);         /* one less element, new total */
    if (alpha <= 1) return;}  /* difference of 0 or 1 is perfect, stop search */

 sum2 = a[0] + a[1];                            /* sum of largest two numbers */
 if (n == 5)                   /* when only 4 numbers are left, KK is optimal */
   {diff2 = sum2 - a[2];                  /* difference of 2 largest elements */
    if (diff2 < a[3]) difference = a[3] - a[4] - diff2;       /* a[3] biggest */
    else difference = diff2 - a[3] - a[4];                /* diff2 is biggest */
    if (difference < 0) difference = - difference;     /* take absolute value */
    if (difference < alpha)            /* better than best previous partition */
      alpha = difference;                      /* reset alpha to better value */
    return;}                                     /* stop searching and return */

 rest = total - a[0] - a[1];          /* sum of all elements except first two */
 if (sum2 >= rest)                       /* if largest element >= sum of rest */
      {if (sum2 - rest < alpha)               /* best better than best so far */
         alpha = sum2 - rest;                  /* reset alpha to better value */
      return;}            /* if difference is larger than rest, so is the sum */

 a[1] = sum2;                    /* sum of two largest is new largest element */
 ckk (a+1, n-1, total);             /* call on subarray with one less element */
 a[1] = sum2 - a[0];                       /* restore array to previous state */

 return;}

/* declares function type  */

long long part (long long [MAX_RNP2009], int, long long, int, int);

  /* INEX is the recursive subroutine that actually searches an
     inclusion-exclusion binary tree to partition a set of numbers k1+k2
     ways. It returns the maximum subset sum in the best partition found so far,
     or minmax if it doesn't find a better partition. */

long long inex
(long long num[MAX_RNP2009],             /* complete set of numbers to be partitioned */
 int n,                                 /* number of numbers to be partitioned */
 long long sum,                    /* sum of all the numbers being partitioned */
 int k,                   /* number of ways set will eventually be partitioned */
 int next,     /* index of next number in NUM array to be included or excluded */
 long long x[MAX_RNP2009],                                 /* numbers included so far */
 int nx,                                          /* number of number included */
 long long y[MAX_RNP2009],                                 /* numbers excluded so far */
 int ny,                                         /* number of numbers excluded */
 long long subsum,                           /* sum of numbers included so far */
 long long rest,                /* sum of numbers not yet included or excluded */
 long long lower,                        /* lower bound on sum of first subset */
 long long upper,                        /* upper bound on sum of first subset */
 long long perfect,                 /* maximum subset sum in perfect partition */
 long long best,      /* maximum subset sum of best k1+k2-way partition so far */
 int top)   /* 1 if this is the top-level partition, 0 if it is a subpartition */

{
  long long newsum;           /* subset sum with an additional number included */
  long long sol1, sol2;          /* maximum subset sums from each subpartition */
  long long newsol;                                /* maximum of sol1 and sol2 */
  long long result;           /* maximum subset sum returned from subpartition */

  if (next == n)    /* set is complete, no more elements to include or exclude */
    {if (k % 2 == 1)                 /* partitioning set an odd number of ways */
        {sol2 = part (y, ny, sum - subsum, k-1, 0); /* optimally partition remainder */
          if (sol2 <= subsum) newsol = subsum;     /* largest of 3 subset sums */
          else newsol = sol2;                  /* largest of three subset sums */
          if (top == 1 && newsol < minmax)           /* better global solution */
            minmax = newsol;                    /* revise best global solution */
          return (newsol);}                       /* return maximum subset sum */
      else                                         /* partitioning set in half */
        {if (nx <= ny)             /* process subset with fewer elements first */
            {sol1 = part (x, nx, subsum, k/2, 0);    /* partition first subset */
              if (sol1 < minmax)             /* could lead to better partition */
                {sol2 = part (y, ny, sum - subsum, k/2, 0); /* partition complement */
                  if (sol1 <= sol2) newsol = sol2;   /* find larger subset sum */
                  else newsol = sol1;     /* largest subset sum of 2 solutions */
                  if (top == 1 && newsol < minmax)   /* better global solution */
                    minmax = newsol;      /* update best solution found so far */
                  return (newsol);}               /* return maximum subset sum */
              else return (minmax);}            /* didn't find better solution */

          else                       /* k1 = k2 and nx > ny, partition y first */
            {sol2 = part (y, ny, sum - subsum, k/2, 0); /* partition complement*/
              if (sol2 < minmax)    /* could possibly lead to better partition */
                {sol1 = part (x, nx, subsum, k/2, 0);       /* part 1st subset */
                  if (sol1 <= sol2) newsol = sol2;   /* find larger subset sum */
                  else newsol = sol1;                /* find larger subset sum */
                  if (top == 1 && newsol < minmax)   /* better global solution */
                    minmax = newsol;           /* update cost of best solution */
                  return (newsol);}         /* return maximum subset sum found */
              else return (minmax);}}}             /* no better solution found */

  newsum = subsum + num[next];            /* subset sum including next element */
  if (newsum <= upper)        /* including next element is still within bounds */
    {x[nx] = num[next];                       /* include next number in subset */
      result = inex (num, n, sum, k, next + 1, x, nx+1, y, ny, newsum, rest-num[next],
                     lower, upper, perfect, best, top);
      if (result < best)             /* found solution better than best so far */
        {if (result == perfect) return (result); /*found optimal solution, quit*/
          best = result;                   /* reset best solution found so far */
          if (k % 2 == 1) lower = sum - (k - 1) * (minmax - 1);
          else lower = sum - k / 2 * (minmax - 1);}}

  if (subsum + rest - num[next] >= lower) /* can reach lower bound excluding next number */
    {y[ny] = num[next];                     /* exclude next number from subset */
      result = inex (num, n, sum, k, next + 1, x, nx, y, ny+1, subsum, rest-num[next],
                   lower, upper, perfect, best, top);
      if (result < best) best = result;}
  return (best);}

/* PART takes an array NUM of integers sorted in decreasing order, the length N
   of the vector, the sum of all the numbers, and a number K of subsets to
   partition it into.  It computes the optimal k-way partition of the numbers if
   all subset sums are less than MINMAX, returning either this value, or MINMAX
   if no better solution was found. It doesn't do any real work, but just calls
   either CKK or INEX to do the actual partitions, depending on whether the
   partition is two ways or more ways. If the last argument, TOP, is equal to
   one, then this function is allowed to reset MINMAX if a better solution is
   found. */

long long part (long long num[MAX_RNP2009], int n, long long sum, int k, int top)

{long long lower;                       /* lower bound of top-level subset sum */
  long long upper;                      /* upper bound of top-level subset sum */
  long long perfect;                /* maximum subset sum in perfect partition */
  long long best;      /* largest subset sum in best partition of argument set */
  int i;                                                      /* utility index */
  long long x[MAX_RNP2009];                                /* numbers in first subset */
  long long y[MAX_RNP2009];                  /* numbers in complement of first subset */


  nodes[k]++;                       /* count number of partitions of each type */

  if (n <= k) return (num[0]);    /* no more sets than numbers, return largest */
  best = kk (num, n, k, sum);    /* initially, Karmarkar-Karp solution is best */
  if (n <= k + 2) return (best);                /* KK is optimal if k <= n + 2 */

  if (k == 2)                                             /* two-way partition */
    {for (i = 0; i < n; i++) /* copy array into temporary since CKK modifies it*/
        A[i] = num[i];
      alpha = sum;               /* set initial difference to very large value */
      ckk (A, n, sum);                    /* compute optimal two-way partition */
      return ((sum + alpha) / 2);}                /* return maximum subset sum */

  if (sum % k == 0) perfect = sum / k;                /* largest subset sum of */
  else perfect = sum / k + 1;                             /* perfect partition */
  if (num[n-1] > perfect) perfect = num[n-1];     /* can't be < largest number */
  if (best == perfect) return (best);      /* kk solution is perfect partition */
  if (top == 1) minmax = best;              /* top-level call, set best so far */

  if (k % 2 == 0)       /* set is to be divided into an even number of subsets */
    {upper = sum / 2;       /* upper bound is half the total sum, rounded down */
      lower = sum - k/2 * (minmax - 1);}        /* lower bound on first subset */
  else                   /* set is to be divided into an odd number of subsets */
    {upper = sum / k;                       /* upper bound on first subset sum */
      lower = sum - (k-1) * (minmax - 1);}      /* lower bound on first subset */

  best = inex (num,n,sum,k,0,x,0,y,0,(long long)0,sum,lower,upper,perfect,best,top);
                 /* search inclusion-exclusion binary tree for better solution */
  if (top == 1) return (minmax); /* top-level call, return global best solution found */
  else return (best);}           /* otherwise return local best solution found */


uint64_t executeRNP2009(const partition::PartitionProblem &problem, ProblemStats &stats) {

  const int &K = problem.K;
  if (K > MAXK) {
    printf("K out of range, MAXK=%d\n", MAXK);
    exit(0);
  }
  N = problem.N;

  memcpy(sorted, problem.S, N * sizeof(int64_t));
  sumall = problem.sum;

  sortdown (sorted, N);                  /* sort numbers in increasing order */

  minmax = sumall; /* initialize best solution to ridiculously large value */

  minmax = part (sorted, N, sumall, K, 1);     /* optimally partition K ways */

  //printf("%ld\n",minmax);
  return minmax;
  //printf ("%3d %12ld\n", trial, minmax);

}

} // end namespace



