 #include <stdio.h>                                    /* standard I/O library */
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
using std::cout;
using std::endl;

#include "RNP.hpp"


namespace partition {


#define MAXK 12         /* maximum number of subsets to partition values into */
#define MAX_IRNP 60                 /* maximum number of values being partitioned */
#define MAXSETS 32768 /* maximum number of values being partitioned 2^(MAXN/4) */
#define MAXLIST 50000000                            /* maximum size of cdlists */
#define INITSEED  13070                    /* initial random seed, in decimal */
#define A 25214903917                                           /* multiplier */
#define C 11                                             /* additive constant */
#define MASK 281474976710655                                      /* 2^{48}-1 */

int K;                            /* number of sets to partition numbers into */
int N;                      /* problem size: number of values in number array */
int TRIALS;                                        /* number of trials to run */
int64_t seed;                                 /* current random number seed */
int64_t MAXNUM;                                    /* maximum integer value */

int64_t sorted[MAX_IRNP];        /* original numbers sorted in decreasing order */
int64_t sumall;                              /* sum of all original numbers */
int64_t minmax;         /* largest subset sum in best solution found so far */

int64_t ckkcalls;            /* number of calls to ckk for 2-way partitions */
int64_t cga3calls;          /* number of calls to cga3 for 3-way partitions */
int64_t cga4calls;          /* number of calls to cga4 for 4-way partitions */
int64_t cga5calls;          /* number of calls to cga4 for 4-way partitions */
int64_t nodes2;                          /* number of SS two-way partitions */
int64_t nodes3;                           /* number of three-way partitions */
int64_t nodes4;                            /* number of four-way partitions */
int64_t nodes5;                            /* number of five-way partitions */

void cga3 (int next, int64_t big, int64_t med, int64_t small);
void cga4 (int next, int64_t huge, int64_t big, int64_t med, int64_t small);
void cga5 (int next, int64_t huge, int64_t big, int64_t med, int64_t small, int64_t tiny);
typedef struct element {           /* structure consists of a set and its sum */
  int64_t sum;
  int64_t set;
} SetArray;

SetArray tsets[MAXSETS];         /* temporary array of sets for merge sorting */
int64_t tsums[MAXSETS];        /* temporary array of sums for merge sorting */

int nextsum;              /* index of next empty position in subset sum array */

typedef struct heapelement{  /* heap of sums from a and b in increasing order */
  int x;                                         /* index of element in asums */
  int y;                                         /* index of element in bsums */
  int64_t sum;}                                      /* sum of the elements */
 Heap;

int64_t numbers[MAX_IRNP];  /* global array of numbers to be partitioned by CKK */

/* global data structures used by TWOPART */

int64_t asums[MAXSETS];              /* subset sums from numbers in A array */
int64_t bsums[MAXSETS];              /* subset sums from numbers in B array */
int64_t csums[MAXSETS];              /* subset sums from numbers in C array */
int64_t dsums[MAXSETS];              /* subset sums from numbers in D array */
int nasums, nbsums, ncsums, ndsums; /* number of sums in each of above arrays */
Heap abheap[MAXSETS], cdheap[MAXSETS]; /*heap of subset sums from pairs of sums*/
int abheapsize, cdheapsize;                       /* size of respective heaps */
int64_t alpha;          /* best difference found for 2-way CKK partitioning */

/* global data structures used by CGA3, CGA4, and CGA5 */

int64_t e[MAX_IRNP];                               /* numbers to be partitioned */
int64_t beta;                                        /* best solution found */
int N345;                    /* number of elements to be partitioned with CGA */

/* INIT initializes the it's argument array with random numbers from zero to
   2^{31}-1.  It returns their sum. */

int64_t init (int64_t a[MAX_IRNP], int n)

{int i;                                                   /* index into array */
 int64_t sum;                                         /* sum of all numbers */

 sum = 0;                                    /* initialize sum of all numbers */
 for (i = 0; i < n; i++)                         /* for each element of array */
   {seed = (A * seed + C) & MASK;             /* next seed in random sequence */

     a[i] = seed;                       /* random value from zero to 2^{48}-1 */
     sum += a[i];}                              /* compute sum of all numbers */
 return (sum);}                                      /* return sum of numbers */

/* SORTINTS takes an array of integers, and the length of the array, and sorts
   the array in increasing order using selection sort. */

void sortints (int64_t a[MAX_IRNP], int n)

{int i,j;                                   /* indices into array for sorting */
 int64_t temp;                              /* temporary value for swapping */

 for (i = 0; i < n-1; i++)                        /* pointer to first element */
   for (j = i+1; j < n; j++)                     /* pointer to second element */
     if (a[i] > a[j])                            /* elements are out of order */
       {temp = a[i];
        a[i] = a[j];                                             /* swap them */
        a[j] = temp;}}

/* SORTLONGS takes an array of longs, and the length of the array, and sorts the
   array in decreasing order using selection sort. */

void sortlongs (int64_t a[MAXK], int k)

{int i,j;                                   /* indices into array for sorting */
 int64_t temp;                              /* temporary value for swapping */

 for (i = 0; i < k-1; i++)                        /* pointer to first element */
   for (j = i+1; j < k; j++)                     /* pointer to second element */
     if (a[i] < a[j])                            /* elements are out of order */
       {temp = a[i];
        a[i] = a[j];                                             /* swap them */
        a[j] = temp;}}

/* INSERT takes an array index FIRST, a new subset sum VECTOR, and an array A of
   vectors sorted in decreasing order of largest elements, and modifies array A
   by inserting the new vector in sorted order by largest element. */

void insert (int first, int n, int k, int64_t vector[MAXK], int64_t vects[MAX_IRNP][MAXK])

{int i;                                  /* index into array of subpartitions */
 int j;                                        /* index to elements of vector */

 for (i = first; i < n-1; i++)         /* insert new partition in sorted list */
   if (vector[0] < vects[i+1][0])          /* haven't found correct place yet */
     for (j = 0; j < k; j++)                    /* for each element of vector */
       vects[i][j] = vects[i+1][j];     /* copy current partition up in order */
   else break;                              /* found correct place, exit loop */

 for (j = 0; j < k; j++)                        /* for each element of vector */
   vects[i][j] = vector[j];}              /* insert new subpartition in order */

/* KK takes an array A of integers, sorted in increasing order, its length N,
   the number of partitions K, and the sum of all the numbers, and returns the
   maximum subset sum in the KK approximation of the best K-way partition. */

int64_t kk (int64_t a[MAX_IRNP], int n, int k, int64_t sum)

{int64_t vects[MAX_IRNP][MAXK];           /* array of vectors N long and K wide */
 int64_t vector[MAXK];                  /* combination of first two vectors */
 int i, j;                                                 /* utility indices */
 int64_t max;                       /* absolute value of maximum subset sum */

 if (k >= n) return ((int64_t) a[n-1]); /*largest subset sum is largest number*/

 for (i = 0; i < n; i++)                   /* create initial array of vectors */
   {vects[i][0] = a[n-i-1];              /* copy to vects in decreasing order */
   for (j = 1; j < k; j++)
     vects[i][j] = 0;}

 for (i = 0; i < n-1; i++)
   {for (j = 0; j < k; j++)                            /* for each subset sum */
     vector[j] = vects[i][j] + vects[i+1][k-j-1];/* combine particular subsets*/
   sortlongs (vector, k);              /* sort the vector in decreasing order */
   for (j = 0; j < k; j++)                             /* for each subset sum */
     vector[j] = vector[j] - vector[k-1];     /* subtract smallest subset sum */
   insert (i+1, n, k, vector, vects);}   /* insert new vector in sorted order */

 for (j = 0; j < k; j++)       /* compute largest subset sum from differences */
   sum = sum - vector[j];                    /* subtract relative subset sums */
 max = (sum / k) + vector[0];                  /* absolute maximum subset sum */
 return (max);}

/* GENSETS takes an array of integers, an array of subsets to fill, the first
   index and last indices into the integer array, the sum of the numbers
   included so far, and the characteristic function of the elements included so
   far. It generates all subsets of the remaining numbers in the NUMbers array
   from index FIRST to index LAST, storing them in the SET array.  It places the
   numbers in consecutive locations indexed by the global variable NEXTSUM, leaving
   in NEXTSUM the number of sets generated. */

void gensets (int64_t num[MAX_IRNP],
              SetArray set[MAXSETS],
              int first, int last,
              int64_t cursum, int64_t curset)

{if (first > last)                                        /* set is completed */
  {set[nextsum].sum = cursum;     /* store subset sum in array of subset sums */
  set[nextsum++].set = curset;}       /* store characteristic function of set */
 else                                                /* set not yet completed */
   {gensets (num, set, first+1, last, cursum, curset);     /* exclude next element */
   gensets (num, set, first+1, last, cursum + num[first], curset | ((int64_t) 1 << first));}}

/* SORTSETS sorts an array of subset sums, and the corresponding subsets
   themselves, in increasing order using mergesort.  It takes the indices of the
   FIRST and LAST elements to be sorted, and sorts those elements in increasing
   order.  It uses a temporary array to do the merge step.  */

void sortsets (SetArray sets[MAXSETS], int first, int last)

{int middle;                                      /* index to middle of array */
 int head1, head2; /* indices of current heads of sorted lists for merge step */
 int tail;                       /* index of next available position in tempa */
 SetArray temp;                 /* temporary pair of set amd sum for swapping */

 if (last - first <= 1 )                 /* pointers are the same or adjacent */
   {if (sets[first].sum > sets[last].sum)     /* two subsets are out of order */
     {temp = sets[first];
     sets[first] = sets[last];                                /* swap subsets */
     sets[last] = temp;}
   return;}                                         /* subarray is now sorted */

 middle = (first + last) / 2;             /* compute middle point of subarray */

 sortsets (sets, first, middle);         /* recursively merge sort first half */
 sortsets (sets, middle+1, last);       /* recursively merge sort second half */

 head1 = first;                                         /* head of first list */
 head2 = middle+1;                                     /* head of second list */
 tail = first;               /* first empty position in sorted temporary list */

 while (head1 <= middle && head2 <= last)    /* neither of the lists is empty */
   if (sets[head1].sum <= sets[head2].sum)   /* head of first list is smaller */
     tsets[tail++] = sets[head1++]; /* grab smaller element and advance head1 */
   else
     tsets[tail++] = sets[head2++]; /* grab smaller element and advance head2 */

 if (head2 > last)                         /* second list was exhausted first */
   while (head1 <= middle)                   /* until first list is exhausted */
     tsets[tail++] = sets[head1++];        /* copy first list to end of merge */

 head1 = first;                        /* reset pointer to head of both lists */
 while (head1 < head2) {                           /* until end of second list */
   sets[head1] = tsets[head1]; head1++;}} /* copy merged list back into original array*/

/* SORTSUMS sorts an array of subset sums in increasing order using mergesort.
   It takes the indices of the FIRST and LAST elements to be sorted, and sorts
   those elements in increasing order.  It uses a temporary array to do the
   merge step.  */

void sortsums (int64_t sums[MAXSETS], int first, int last)

{int middle;                                      /* index to middle of array */
 int head1, head2; /* indices of current heads of sorted lists for merge step */
 int tail;                       /* index of next available position in tempa */
 int64_t temp;                /* temporary pair of set amd sum for swapping */

 if (last - first <= 1 )                 /* pointers are the same or adjacent */
   {if (sums[first] > sums[last])     /* two subsets are out of order */
     {temp = sums[first];
     sums[first] = sums[last];                                /* swap subsets */
     sums[last] = temp;}
   return;}                                         /* subarray is now sorted */

 middle = (first + last) / 2;             /* compute middle point of subarray */

 sortsums (sums, first, middle);         /* recursively merge sort first half */
 sortsums (sums, middle+1, last);       /* recursively merge sort second half */

 head1 = first;                                         /* head of first list */
 head2 = middle+1;                                     /* head of second list */
 tail = first;               /* first empty position in sorted temporary list */

 while (head1 <= middle && head2 <= last)    /* neither of the lists is empty */
   if (sums[head1] <= sums[head2]) {           /* head of first list is smaller */
     tsums[tail] = sums[head1]; /* grab smaller element and advance head1 */
     tail++;
     head1++;
   } else {
     tsums[tail] = sums[head2]; /* grab smaller element and advance head2 */
     tail++;
     head2++;
   }
 if (head2 > last)                         /* second list was exhausted first */
   while (head1 <= middle)                   /* until first list is exhausted */
     tsums[tail++] = sums[head1++];        /* copy first list to end of merge */

 head1 = first;                        /* reset pointer to head of both lists */
 while (head1 < head2) {                            /* until end of second list */
   sums[head1] = tsums[head1];head1++;}} /* copy merged list back into original array*/

/* GENSUMS takes an array X of numbers, the index of the LAST number, and
   generates an array SUMS of all sums that can be achieved by adding together
   numbers from X.  Each number in X can only be used at most once.  To allow
   the recursion, additional arguments include an index NEXTX to the next
   element of the array, and the SUMSOFAR achieved down this path.  There is
   also a global pointer NEXTSUM to the next empty element of the SUMS array. At
   the end, it is equal to the number of sums created. */

void gensums (
		int64_t nums[MAX_IRNP],                          /* array of original numbers */
int64_t sums[MAXSETS],                                   /* array of sums */
int nextx,                            /* pointer to next element of X array */
int64_t sumsofar,             /* sum so far of elements in current subset */
int last                                 /* index of last element in array */
)
{
	if (nextx == last)                          /* reached last element of array */
   {sums[nextsum++] = sumsofar;                     /* don't add last element */
    sums[nextsum++] = sumsofar + nums[nextx];}        /* add last element to sum */
  else {gensums (nums, sums, nextx+1, sumsofar, last);     /* don't add last element */
    gensums (nums, sums, nextx+1, sumsofar + nums[nextx], last);}} /* add last element */

/* CKK takes an array of int64_t numbers A, its length N, and their TOTAL sum,
   and finds the best partition, leaving the resulting difference in the global
   variable ALPHA. It runs branch-and-bound, starting with the Karmarkar-Karp
   solution. At each point, the largest two remaining numbers are selected, and
   replaced with either their difference or their sum, representing assigning
   them to different sets or the same set, respectively. */

void ckk (
int64_t a[MAX_IRNP],                                        /* array of numbers */
int n,                                         /* number of elements in array */
int64_t total)                                        /* sum of all numbers */

{int64_t diff2;                          /* difference of largest 2 numbers */
 int64_t sum2;                                  /* sum of largest 2 numbers */
 int64_t difference;                             /* difference of partition */
 int64_t b[MAX_IRNP];                   /* new copy of list for recursive calls */
 int64_t rest;                        /* sum of all elements except first 2 */
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
 ckk (a+1, n-1, total);          /* call on subarray with one less element */
 a[1] = sum2 - a[0];                       /* restore array to previous state */

 return;}

/* TWOPART takes an array of NUMberS, the Number of numbers in the array, and
   their SUM, and optimally partitions them two ways using the Schroeppel and
   Shamir algorithm.  It returns the larger subset sum in the optimal
   partition. */

int64_t twopart (int64_t nums[MAX_IRNP], int n, int64_t sum)

{int i;                                                      /* utility index */
 int64_t best;                           /* best partition difference found */
 int64_t target1, target2;       /* two target values, equal if sum is even */
 int dpointer;                 /* index to largest number in dsums <= target2 */
 int child, parent;              /* indices of parent and child nodes in heap */
 Heap next;                                     /* new element to add to heap */
 int64_t newsum;                         /* new subset sum being considered */
 int64_t complement;                                        /* sum - newsum */
 int64_t newmax;                        /* maximum of newsum and complement */

 if (n < 17)               /* small number of elements, run CKK instead of SS */
   {alpha = MASK;               /* initial value of best difference, 2^{48}-1 */
   ckkcalls++;                                /* count number of calls to ckk */
   for (i = 0; i < n; i++) /* */
     numbers[i] = nums[n-1-i]; /* copy numbers in decreasing order into global array*/
   ckk (numbers, n, sum);        /* run CKK to find best 2-way partition difference */
   return ((sum + alpha) / 2);}  /* return larger subset of optimal partition */

 target1 = sum / 2;                                   /* smaller target value */
 target2 = sum - target1;                              /* larger target value */

 nodes2++;                           /* count number of SS two-way partitions */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, asums, 0, (int64_t) 0, n/4-1);        /* generate sums of elements */
 nasums = nextsum;                                /* number of sums generated */
 sortsums (asums, 0, nasums-1);            /* merge sort array of sums from a */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, bsums, n/4, (int64_t) 0, n/2-1);      /* generate sums of elements */
 nbsums = nextsum;                                /* number of sums generated */
 sortsums (bsums, 0, nbsums-1);            /* merge sort array of sums from b */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, csums, n/2, (int64_t) 0, n*3/4-1);    /* generate sums of elements */
 ncsums = nextsum;                                /* number of sums generated */
 sortsums (csums, 0, ncsums-1);            /* merge sort array of sums from c */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, dsums, n*3/4, (int64_t) 0, n-1);      /* generate sums of elements */
 ndsums = nextsum;                                /* number of sums generated */
 sortsums (dsums, 0, ndsums-1);            /* merge sort array of sums from d */

 for (i = 0; i < nasums; i++)                                 /* make ab heap */
   {abheap[i].x = i;                         /* aptr pointer is same as index */
    abheap[i].y = 0;              /* smallest value of bsums is first element */
    abheap[i].sum = asums[i];        /* asums[x] + bsums[0], but bsums[0] = 0 */
    if (abheap[i].sum > target2) break;}     /* sum too big, done making heap */
 abheapsize = i;                      /* number of elements in heap initially */

                /* make CD maxheap, only with elements <= larger target value */
 dpointer = ndsums - 1;                    /* index to largest value in dsums */
 cdheapsize = 0;                                   /* no elements in heap yet */
 for (i = 0; i < ncsums; i++)            /* process csums in increasing order */
   {if (csums[i] > target2) break;          /* no more entries to add to heap */
   else cdheapsize++;                 /* index of next empty location in heap */
   next.x = i;                                           /* cpointer is index */
   while (csums[i] + dsums[dpointer] > target2) dpointer--;   /* sums too big */
   next.y = dpointer;                     /* largest value of dsums <= target */
   next.sum = csums[i] + dsums[dpointer];          /* sum of new heap element */
   child = cdheapsize-1;                        /* index of last node in heap */
   while (child > 0)                   /* bubble child up to correct location */
     {parent = (child - 1) / 2;           /* index in heap of parent of child */
       if (next.sum > cdheap[parent].sum)   /* new element is too low in heap */
         {cdheap[child] = cdheap[parent];                 /* move parent down */
           child = parent;}                       /* parent is new child node */
       else break;}             /* found correct location of new heap element */
   cdheap[child] = next;}                              /* insert it into heap */

 best = minmax;              /* only interested in solutions better than this */

 while (abheapsize > 0 && cdheapsize > 0)  /* there are still valid sums left */
   {newsum = abheap[0].sum + cdheap[0].sum; /* next ab number + next cd number*/
   if (newsum == target1 || newsum == target2)     /* found perfect partition */
     return (target2);                            /* return larger subset sum */
   complement = sum - newsum;              /* subset sum of complementary set */
   if (complement <= newsum) newmax = newsum;     /* newsum is larger of pair */
   else newmax = complement;                          /* complement is larger */
   if (newmax < best) best = newmax;  /* new solution better than best so far */

   if (newsum < target1)                   /* construct next element for heap */
     {if (abheap[0].y + 1 < nbsums)             /* still more b pointers left */
         {next.x = abheap[0].x;                   /* x pointer stays the same */
           next.y = abheap[0].y + 1;          /* set to next larger b pointer */
           next.sum = asums[next.x] + bsums[next.y]; /*sum of new heap element*/
           if (next.sum > target2)          /* this value is already too large */
             {next = abheap[abheapsize-1]; /*replace top with last element of heap*/
               abheapsize--;}}                  /* heap has one fewer element */
       else                                  /* remove this element from heap */
         {abheapsize--;                         /* heap has one fewer element */
           next = abheap[abheapsize];} /*replace top with last element of heap*/

                       /* percolate new top of heap down to its correct place */
       parent = 0;
       while (2 * parent + 1 < abheapsize)   /* still at least one child left */
         {if (2 * parent + 2 == abheapsize)        /* node only has one child */
             child = 2 * parent + 1;
           else if (abheap[2*parent+1].sum < abheap[2*parent+2].sum) child = 2*parent+1;
           else child = 2 * parent + 2;                 /* find smaller child */
           if (next.sum <= abheap[child].sum) break;           /* found place */
           else
             {abheap[parent] = abheap[child];    /* bubble child up to parent */
               parent = child;}}                   /* child is now new parent */
       abheap[parent] = next;}      /* insert next element in parent position */

   else                  /* sum is greater than target, remove top of cd heap */
     {if (cdheap[0].y-1 >= 0)                   /* still more d pointers left */
         {next.y = cdheap[0].y - 1;          /* set to next smaller d pointer */
           next.x = cdheap[0].x;                     /* x pointer is the same */
           next.sum = csums[next.x] + dsums[next.y];}
       else                                  /* remove this element from heap */
         {cdheapsize--;                         /* heap has one fewer element */
           next = cdheap[cdheapsize];} /*make new element last element of heap*/
                       /* percolate new top of heap down to its correct place */
       parent = 0;
       while (2 * parent + 1 < cdheapsize)   /* still at least one child left */
         {if (2 * parent + 2 == cdheapsize)        /* node only has one child */
             child = 2 * parent + 1;
           else if (cdheap[2*parent+1].sum > cdheap[2*parent+2].sum) child = 2*parent+1;
           else child = 2 * parent + 2;                  /* find larger child */
           if (next.sum >= cdheap[child].sum) break;           /* found place */
           else
             {cdheap[parent] = cdheap[child];
               parent = child;}}                   /* child is now new parent */
       cdheap[parent] = next;}}     /* insert next element in parent position */
 return (best);}

/* PART takes an array NUM of integers sorted in increasing order, the length N
   of the vector, the sum of all the numbers, and a number K of subsets to
   partition it into.  It computes the optimal k-way partition of the numbers if
   all subset sums are less than MINMAX, returning either this value, or MINMAX
   if no better solution was found. If the last argument, TOP, is equal to one,
   then this function is allowed to reset MINMAX if a better solution is
   found. */

int64_t part (int64_t num[MAX_IRNP], int n, int64_t sum, int k, int top)

{int64_t lower;                       /* lower bound of top-level subset sum */
  int64_t upper;                      /* upper bound of top-level subset sum */
  int64_t perfect;                /* maximum subset sum in perfect partition */

  SetArray a[MAXSETS];    /* array of large sets from first quarter of numbers */
  SetArray b[MAXSETS];   /* array of large sets from second quarter of numbers */
  SetArray c[MAXSETS];    /* array of large sets from third quarter of numbers */
  SetArray d[MAXSETS];    /* array of large sets from third quarter of numbers */
  int nasets, nbsets, ncsets, ndsets;          /* number of sets in each array */

  Heap abheap[MAXSETS];       /* heap of sets from combination of A and B sets */
  Heap cdheap[MAXSETS];       /* heap of sets from combination of C and D sets */
  int abheapsize, cdheapsize;               /* number of elements in each heap */
  int dpointer; /* pointer to current element in d array, when building CD heap*/
  Heap cdlist[MAXLIST];           /* list of sets already removed from CD heap */
  int cdlistsize;             /* number of subsets in list popped from CD Heap */
  Heap newElement;                                 /* new element being added to heap */

  int64_t x[MAX_IRNP], y[MAX_IRNP]; /* arrays to contain first subset and complement */
  int nx, ny;                     /* number of numbers in each of these arrays */
  int first;                                /* index of first subset in cdlist */
  int firstminusone;               /* index just before first subset in cdlist */
  int next;                         /* index of first empty location in cdlist */
  int64_t topsum;        /* sum of combination subset from top of both heaps */
  int index;                                              /* index into cdlist */
  int64_t sol1, sol2; /*maximum subset sums from subpartitions of two subsets*/
  int64_t newsol;                                /* maximum of sol1 and sol2 */
  int64_t set;           /* characteristic function of set being constructed */
  int64_t subsum;                              /* subset sum of first subset */
  int64_t compsum;               /* subset sum of complement of first subset */
  int k1, k2; /*number of subsets to partition first subset and complement into*/
  int parent, child;              /* indices of parent and child nodes in heap */
  int i;                                                      /* utility index */
  int64_t best;      /* largest subset sum in best partition of argument set */

  if (k == 3 && n < 13)             /* special case small three-way partitions */
      {beta = minmax; /* initial upper bound is max sum in best solution so far*/
      cga3calls++;                            /* count number of calls to cga3 */
      for (i = 0; i < n; i++)
        e[i] = num[n-1-i]; /*copy numbers in decreasing order into global array*/
      N345 = n;  /* global number of elements to be partitioned for efficiency */
      cga3 (1, (int64_t) e[0], (int64_t) 0, (int64_t) 0); /*find solution*/
      return (beta);}          /* return largest subset in best solution found */

  if (k == 4 && n < 15)              /* special case small four-way partitions */
    {beta = minmax;      /* initial upper bound is max of best solution so far */
      cga4calls++;                            /* count number of calls to cga4 */
      for (i = 0; i < n; i++)
        e[i] = num[n-1-i]; /*copy numbers in decreasing order into global array*/
      N345= n;   /* global number of elements to be partitioned for efficiency */
      cga4 (1, (int64_t) e[0], (int64_t) 0, (int64_t) 0, (int64_t) 0);
      return (beta);}          /* return largest subset in best solution found */

  if (k == 5 && n < 17)              /* special case small five-way partitions */
    {beta = minmax;      /* initial upper bound is max of best solution so far */
      cga5calls++;                            /* count number of calls to cga5 */
      for (i = 0; i < n; i++)                      /* for each original number */
        e[i] = num[n-1-i];  /* copy them in decreasing order into global array */
      N345 = n;  /* global number of elements to be partitioned for efficiency */
      cga5(1, (int64_t) e[0], (int64_t) 0, (int64_t) 0, (int64_t) 0, (int64_t) 0);
      return (beta);}          /* return largest subset in best solution found */

  best = kk (num, n, k, sum);    /* initially, Karmarkar-Karp solution is best */
  if (n <= k + 2) return (best);         /* kk is always optimal if k <= n + 2 */
  if (sum % k == 0) perfect = sum / k;                /* largest subset sum of */
  else perfect = sum / k + 1;                             /* perfect partition */
  if (num[n-1] > perfect) perfect = num[n-1];     /* can't be < largest number */
  if (best == perfect) return (best);      /* kk solution is perfect partition */
  if (top == 1) minmax = best;              /* top-level call, set best so far */

  if (k == 3) nodes3++;                /* count number of three-way partitions */
  else if (k == 4) nodes4++;            /* count number of four-way partitions */
  else if (k == 5) nodes5++;            /* count number of four-way partitions */

  k1 = k / 2;       /* number of subsets first subset will be partitioned into */
  k2 = k - k1; /* number of subsets remaining elements will be partitioned into*/
  upper = k1 * (minmax - 1);                /* upper bound on first subset sum */
  if (sum % k <= k1) lower = k1 * (sum / k) + sum % k; /* lower bound based on */
  else lower = k1 * (sum / k + 1);                        /* perfect partition */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets (num, a, 0, n/4-1, 0, 0);  /* combinations of 1st quarter of numbers */
  nasets = nextsum;                             /* number of subsets generated */
  sortsets (a, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets (num, b, n/4, n/2-1, 0, 0); /*combinations of 2nd quarter of numbers */
  nbsets = nextsum;                             /* number of subsets generated */
  sortsets (b, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets (num, c, n/2, n*3/4-1, 0, 0);/*combinations of 3rd quarter of numbers*/
  ncsets = nextsum;                             /* number of subsets generated */
  sortsets (c, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets (num, d, n*3/4, n-1, 0, 0); /*combinations of 4th quarter of numbers */
  ndsets = nextsum;                             /* number of subsets generated */
  sortsets (d, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  /*construct min heap of initial combinations of sets from array A and array B*/
  for (index = 0; index < nasets; index++)          /* for each set in A array */
    {abheap[index].x = index;                 /* aptr pointer is same as index */
      abheap[index].y = 0;         /* smallest value of bsums is first element */
      abheap[index].sum = a[index].sum; /*a[index].sum + b[0].sum, b[0].sum = 0*/
      if (abheap[index].sum > upper) break;}  /* sum too big, done making heap */
  abheapsize = index;                          /* number of elements in abheap */

  /*construct max heap of initial combinations of sets from array C and array D*/
  dpointer = ndsets - 1;                    /* index to largest value in dsums */
  cdheapsize = 0;                                   /* no elements in heap yet */
  for (index = 0; index < ncsets; index++) /*process csums in increasing order */
    {if (c[index].sum > upper) break;        /* no more entries to add to heap */
      else cdheapsize++;               /* index of next empty location in heap */
      newElement.x = index;                     /* cpointer is index */
      while (c[index].sum + d[dpointer].sum > upper) dpointer--; /* find sums <= upper */
      newElement.y = dpointer;   /* largest value of dsums <= target */
      newElement.sum = c[index].sum + d[dpointer].sum; /* sum of heap element */
      child = cdheapsize-1;                     /* index in heap of child node */
      while (child > 0)                 /* bubble child up to correct location */
        {parent = (child - 1) / 2;         /* index in heap of parent of child */
          if (newElement.sum > cdheap[parent].sum)/*child is too low in heap*/
            {cdheap[child] = cdheap[parent];
              child = parent;}                     /* parent is new child node */
          else break;} /* child is in correct place, add next element to heap */
      cdheap[child] = newElement;}                              /* insert it into heap */

  first = 0;                                         /* no entries in list yet */
  next = 0;                                    /* first empty location in list */
  cdlistsize = 0;                     /* list of cd subsets is initially empty */

  while (abheapsize > 0)                             /* until AB heap is empty */
    /* first remove sets with sums above upper bound from list */
    {while (cdlistsize > 0 && abheap[0].sum + cdlist[first].sum > upper)
        {first = (first + 1) % MAXLIST;             /* remove subset from list */
          cdlistsize--;}                                /* reduce size of list */

      /* add new sets to list with sums within bounds */
      topsum = cdheap[0].sum + abheap[0].sum; /*sum of elements on top of heaps*/
      while (topsum >= lower)                         /* new set within bounds */
        {if (topsum <= upper)                         /* new set within bounds */
            {cdlist[next] = cdheap[0];             /* add heap element to list */
              next = (next + 1) % MAXLIST;     /* add new set to end of ablist */
              cdlistsize++;
              if (cdlistsize > MAXLIST)
                {printf ("error: subset list size exceeded\n"); exit(0);}}

          if (cdheap[0].y > 0)      /* there is another element in this column */
            {newElement.y = cdheap[0].y - 1;    /* replace with next combination in this column */
              newElement.x = cdheap[0].x;
              newElement.sum = c[newElement.x].sum + d[newElement.y].sum;}/*sum of new subset*/
          else                                  /* no more elements in this column */
            {cdheapsize--;                      /* size of heap shrinks by one */
              if (cdheapsize > 0)                   /* CDheap is not yet empty */
                newElement = cdheap[cdheapsize];/*replace top of heap with last member*/
              else break;}         /* abheap is empty, break out of while loop */
          /* percolate new top of heap down to its correct position */
          parent = 0;                    /* initially parent node is root node */
          while (2 * parent + 1 < cdheapsize) /* still at least one child left */
            {if (2 * parent + 2 == cdheapsize) /* node has only one child left */
                child = 2 * parent + 1;
              else if (cdheap[2*parent+1].sum >= cdheap[2*parent+2].sum)
                child = 2 * parent + 1;               /* left child is smaller */
              else child = 2 * parent + 2;           /* right child is smaller */
              if (newElement.sum >= cdheap[child].sum) break; /*found place*/
              else
                {cdheap[parent] = cdheap[child];
                  parent = child;}}                 /* child is now new parent */
          cdheap[parent] = newElement;     /* insert next element in parent position */
          topsum = cdheap[0].sum + abheap[0].sum;}  /* sum of two top elements of heaps */

      if (cdlistsize > 0)/*any subsets left in the list between lower and upper*/
        {if (next == 0) index = MAXLIST - 1; /* index of last element in cdlist*/
          else index = next - 1;               /* process sums in decreasing order */
          if (first == 0) firstminusone = MAXLIST - 1; /* index before first in cdlist */
          else firstminusone = first - 1; /* */
          while (index != firstminusone) /* */
            {subsum = cdlist[index].sum + abheap[0].sum;   /* first subset sum */
              if (subsum >= lower && subsum <= upper)/*subset sum within bounds*/
                {set = c[cdlist[index].x].set |
                       d[cdlist[index].y].set |
                       a[abheap[0].x].set |
                       b[abheap[0].y].set;
                  compsum = sum - subsum;
                  if (k == 3)                   /* partitioning set three ways */
                    {ny = 0;  /* constuct array of numbers not in first subset */
                      for (i = 0; i < n; i++)
                        {if ((set & 1) == 0)  /* number is not in first subset */
                            {y[ny] = num[i];  /* add element to complement set */
                              ny++;}        /* count numbers in complement set */
                          set = set >> 1;} /* shift to get bit for next element*/
                      sol2 = twopart (y, ny, compsum); /* larger of 2-way part */
                      if (sol2 <= subsum) newsol = subsum; /* largest of 3 subset sums */
                      else newsol = sol2;      /* largest of three subset sums */
                      if (newsol < best)    /* better than best local solution */
                        {best = newsol;          /* update best local solution */
                          if (top == 1 && newsol < minmax) /* better global solution */
                            {minmax = newsol;   /* revise best global solution */
                              upper = k1 * (minmax - 1); /* new bound on first subset sum*/
                              if (lower > upper || minmax == perfect) return (best);}}}
                  else                    /* partitioning more than three ways */
                    {ny = nx = 0;   /* constuct arrays of numbers in both sets */
                      for (i = 0; i < n; i++)
                        {if ((set & 1) == 1)  /* number is not in first subset */
                            {x[nx] = num[i];  /* add element to complement set */
                              nx++;}        /* count numbers in complement set */
                          else
                            {y[ny] = num[i];    /* add element to original set */
                              ny++;}        /* count numbers in complement set */
                          set = set >> 1;} /* shift to get bit for next element*/
                      if (nx <= ny)                  /* fewer numbers in x set */
                        {if (k1 == 2) sol1 = twopart (x, nx, subsum); /* 2-way partition */
                          else sol1 = part (x, nx, subsum, k1, 0);/*partition first subset*/
                          if (sol1 < minmax) /* could lead to better partition */
                            {if (k2 == 2) sol2 = twopart (y, ny, compsum);
                              else sol2 = part (y, ny, compsum, k2, 0);
                              if (sol1 <= sol2) newsol = sol2; /* find larger subset sum */
                              else newsol = sol1; /* largest subset sum of 2 solutions */
                              if (newsol < best)  /* better than best local solution */
                                {best = newsol;  /* update best local solution */
                                  if (top == 1 && newsol < minmax)/*better global solution*/
                                    {minmax = newsol;
                                      upper = k1 * (minmax - 1); /*new bound on first set */
                                      if (lower > upper || minmax == perfect)
                                        return (best);}}}} /* found best local solution */

                      else    /* x has more elements than y, partition y k2 ways first */
                        {if (k2 == 2) sol2 = twopart (y, ny, compsum); /* 2-way part */
                          else sol2 = part (y, ny, compsum, k2, 0); /* >2-way part */
                          if (sol2 < minmax)   /* could lead to better partition */
                            {if (k1 == 2) sol1 = twopart (x, nx, subsum);/*part 1st 2 ways*/
                              else sol1 = part (x, nx, subsum, k1, 0); /* part 1st >2 ways*/
                              if (sol1 <= sol2) newsol = sol2; /*find larger subset sum */
                              else newsol = sol1;    /* find larger subset sum */
                              if (newsol < best) /* better than best local solution */
                                {best = newsol;    /* update best local solution */
                                  if (top == 1 && newsol < minmax)/*better global solution*/
                                    {minmax = newsol;   /* update cost of best solution */
                                      upper = k1 * (minmax - 1);/*new bound on first subset*/
                                      if (lower > upper || minmax == perfect)
                                        return (best);}}}}}} /* perfect partition, quit */

              if (index == 0) index = MAXLIST - 1;       /* previous index in cdlist */
              else index = index - 1;}}  /* */

      if (abheap[0].y < nbsets - 1)            /* there is another element in this column */
        {newElement.x = abheap[0].x;          /* replace with next combination in this column */
          newElement.y = abheap[0].y + 1;
          newElement.sum = a[newElement.x].sum + b[newElement.y].sum;} /* sum of new subset */
      else
        {abheapsize--;                           /* decrease size of heap by one */
          if (abheapsize > 0)                           /* ABheap is not yet empty */
            newElement = abheap[abheapsize]; /* replace top of heap with last member*/
          else break;}                                          /* cdheap is empty */
                    /* percolate new top of heap down to its correct position */
      parent = 0;                          /* initially parent node is root node */
      while (2 * parent + 1 < abheapsize)       /* still at least one child left */
        {if (2 * parent + 2 == abheapsize)       /* node has only one child left */
            child = 2 * parent + 1;
          else if (abheap[2*parent+1].sum <= abheap[2*parent+2].sum) /* compare children */
            child = 2 * parent + 1;                        /* left child is larger */
          else child = 2 * parent + 2;                    /* right child is larger */
          if (newElement.sum <= abheap[child].sum) break; /* found right place */
          else
            {abheap[parent] = abheap[child];
              parent = child;}}                          /* child is now new parent */
      abheap[parent] = newElement;}      /* insert next element in parent position */

  return (best);}

/* CGA3 takes the index of the first unassigned number, the sums of the two
 largest subsets, minus the size of the smallest, and the sum of the remaining
 numbers, and finds the best partition for the remaining numbers, leaving the
 resulting difference in beta. */

void cga3 (int next, int64_t big, int64_t med, int64_t small)

{int64_t newNumber;                  /* new number created by adding two together */

 if (next == N345-1)                                       /* one number left */
   {newNumber = small + e[next];        /* place last number in smallest subset sum */
   if (newNumber < beta) {                               /* new best partition found */
     if (newNumber >= big) beta = newNumber;                     /* new number is largest */
     else beta = big;                                 /* big is still largest */
   }
   return;}                                /* solution is complete, backtrack */

 newNumber = small + e[next];     /* place next number in smallest subset sum first */
 if (newNumber < beta)
   {if (newNumber >= big)                              /* new number is now largest */
     cga3 (next+1, newNumber, big, med);         /* recursively search rest of tree */
   else if (newNumber >= med)               /* new number is between big and medium */
     cga3 (next+1, big, newNumber, med);         /* recursively search rest of tree */
   else                                    /* next number is biggest of three */
     cga3 (next+1, big, med, newNumber);         /* recursively search rest of tree */
   if (big >= beta) return;}          /* recursive call found better solution */

 if (med > 0)         /* if both small and med are zero, don't put it in both */
   {newNumber = med + e[next];                  /* put next number in medium subset */
   if (newNumber < beta)
     {if (newNumber >= big)              /* new element plus medium is still medium */
       cga3 (next+1, newNumber, big, small);     /* recursively search rest of tree */
     else                                        /* new subset is now largest */
       cga3 (next+1, big, newNumber, small);     /* recursively search rest of tree */
     if (big >= beta) return;}}       /* recursive call found better solution */

 if (big > 0)              /* if all three are zero, don't place in all three */
   {newNumber = big + e[next];                    /* put next number in largest set */
   if (newNumber < beta) cga3 (next+1, newNumber, med, small);}}   /* search rest of tree */

/* CGA4 takes the index of the first unassigned number, the sums of the two
 largest subsets, minus the size of the smallest, and the sum of the remaining
 numbers, and finds the best partition for the remaining numbers, leaving the
 resulting difference in alpha. */

void cga4 (int next, int64_t huge, int64_t big, int64_t med, int64_t small)

{int64_t newNumber;                  /* new number created by adding two together */

 if (next == N345-1)                                       /* one number left */
   {newNumber = small + e[next];      /* place next number in smallest subset first */
   if (newNumber < beta) {                               /* new best partition found */
     if (newNumber >= huge) beta = newNumber;                    /* new number is largest */
     else beta = huge;                                /* big is still largest */
   }
   return;}                                /* solution is complete, backtrack */

 newNumber = small + e[next]; /* place next number in set with smallest subset first*/
 if (newNumber < beta)
   {if (newNumber >= huge)                             /* new number is now largest */
     cga4 (next+1, newNumber, huge, big, med);   /* recursively search rest of tree */
   else if (newNumber >= big)               /* new number is between big and medium */
     cga4 (next+1, huge, newNumber, big, med);   /* recursively search rest of tree */
   else if (newNumber >= med)               /* new number is between big and medium */
     cga4 (next+1, huge, big, newNumber, med);   /* recursively search rest of tree */
   else                                    /* next number is biggest of three */
     cga4 (next+1, huge, big, med, newNumber);   /* recursively search rest of tree */
   if (huge >= beta) return;}         /* recursive call found better solution */

 if (med > 0)         /* if both small and med are zero, don't put it in both */
   {newNumber = med + e[next];                  /* put next number in medium subset */
   if (newNumber < beta)
     {if (newNumber >= huge)             /* new element plus medium is still medium */
       cga4 (next+1, newNumber, huge, big, small);           /* search rest of tree */
     else if (newNumber >= big)          /* new element plus medium is still medium */
       cga4 (next+1, huge, newNumber, big, small);           /* search rest of tree */
     else cga4 (next+1, huge, big, newNumber, small);        /* search rest of tree */
     if (huge >= beta) return;}}      /* recursive call found better solution */

 if (big > 0)              /* if all three are zero, don't place in all three */
   {newNumber = big + e[next];                    /* put next number in largest set */
   if (newNumber < beta)
     {if (newNumber > huge)                               /* new element is largest */
       cga4 (next+1, newNumber, huge, med, small);           /* search rest of tree */
     else cga4 (next+1, huge, newNumber, med, small);        /* search rest of tree */
     if (huge >= beta) return;}}      /* recursive call found better solution */

 if (huge > 0)
   {newNumber = huge + e[next];                   /* put next number in largest set */
   if (newNumber < beta)
     cga4 (next+1, newNumber, big, med, small);}} /*recursively search rest of tree */

/* CGA5 takes the index of the first unassigned number, the sums of the five
 largest subsets so far, and finds the best partition for the remaining numbers,
 leaving the resulting difference in beta. */

void cga5 (int next, int64_t huge, int64_t big, int64_t med, int64_t small, int64_t tiny)

{int64_t newNumber;                  /* new number created by adding two together */

 if (next == N345 - 1)                                     /* one number left */
   {newNumber = tiny + e[next];       /* place next number in smallest subset first */
   if (newNumber < beta) {                               /* new best partition found */
     if (newNumber >= huge) beta = newNumber;                    /* new number is largest */
     else beta = huge;                                /* big is still largest */
   }
   return;}                                /* solution is complete, backtrack */

 newNumber = tiny + e[next]; /* place next number in set with smallest subset first */
 if (newNumber < beta)
   {if (newNumber >= huge)                             /* new number is now largest */
     cga5 (next+1, newNumber, huge, big, med, small);        /* search rest of tree */
   else if (newNumber >= big)               /* new number is between big and medium */
     cga5 (next+1, huge, newNumber, big, med, small);        /* search rest of tree */
   else if (newNumber >= med)               /* new number is between big and medium */
     cga5 (next+1, huge, big, newNumber, med, small);        /* search rest of tree */
   else if (newNumber >= small)                  /* next number is biggest of three */
     cga5 (next+1, huge, big, med, newNumber, small);        /* search rest of tree */
   else cga5 (next+1, huge, big, med, small, newNumber);     /* search rest of tree */
   if (huge >= beta) return;}         /* recursive call found better solution */

 if (small > 0)       /* if both small and med are zero, don't put it in both */
   {newNumber = small + e[next];                /* put next number in medium subset */
   if (newNumber < beta)
     {if (newNumber >= huge)             /* new element plus medium is still medium */
       cga5 (next+1, newNumber, huge, big, med, tiny);       /* search rest of tree */
     else if (newNumber >= big)          /* new element plus medium is still medium */
       cga5 (next+1, huge, newNumber, big, med, tiny);       /* search rest of tree */
     else if (newNumber >= med)
       cga5 (next+1, huge, big, newNumber, med, tiny);       /* search rest of tree */
     else cga5 (next+1, huge, big, med, newNumber, tiny);    /* search rest of tree */
     if (huge >= beta) return;}}      /* recursive call found better solution */

 if (med > 0)              /* if all three are zero, don't place in all three */
   {newNumber = med + e[next];                    /* put next number in largest set */
   if (newNumber < beta)
     {if (newNumber >= huge)                              /* new element is largest */
       cga5 (next+1, newNumber, huge, big, small, tiny);     /* search rest of tree */
     else if (newNumber >= big)
       cga5 (next+1, huge, newNumber, big, small, tiny);     /* search rest of tree */
     else cga5 (next+1, huge, big, newNumber, small, tiny);  /* search rest of tree */
     if (huge >= beta) return;}}      /* recursive call found better solution */

 if (big > 0)
   {newNumber = big + e[next];                    /* put next number in largest set */
   if (newNumber < beta)
     {if (newNumber >= huge)
       cga5 (next+1, newNumber, huge, med, small, tiny);     /* search rest of tree */
     else cga5 (next+1, huge, newNumber, med, small, tiny);  /* search rest of tree */
     if (huge >= beta) return;}}      /* recursive call found better solution */

 if (huge > 0)
   {newNumber = huge + e[next];
   if (newNumber < beta)
     cga5 (next+1, newNumber, big, med, small, tiny);}}

uint64_t executeRNP(const partition::PartitionProblem &problem, ProblemStats &stats) {

  const int &K = problem.K;
  if (K > MAXK) {
    printf("K out of range, MAXK=%d\n", MAXK);
    exit(0);
  }
  N = problem.N;

  memcpy(sorted, problem.S, N * sizeof(int64_t));
  sumall = problem.sum;

  sortints (sorted, N);                  /* sort numbers in increasing order */

  minmax = sumall; /* initialize best solution to ridiculously large value */

  if (K == 2) {
    cout.flush();
    minmax = twopart(sorted, N, sumall / 2); /*optimal 2-way partition*/
  } else {
    minmax = part(sorted, N, sumall, K, 1); /* optimally partition K ways*/
  }

  //printf("%ld\n",minmax);
  return minmax;
  //printf ("%3d %12ld\n", trial, minmax);

}
}
