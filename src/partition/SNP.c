#include <stdio.h>                                    /* standard I/O library */
#include <math.h>                                     /* mathematical library */
#include <stdlib.h>

#define MAXK 12         /* maximum number of subsets to partition values into */
#define MAXN 64                 /* maximum number of values being partitioned */
#define MAXSETS 65536 /* maximum number of values being partitioned 2^(MAXN/4)*/
#define MAXLIST 100000                             /* maximum size of cdlists */
#define INITSEED  13070                    /* initial random seed, in decimal */
#define ACONST 25214903917                             /* constant multiplier */
#define C 11                                             /* additive constant */
#define MASK 281474976710655                                      /* 2^{48}-1 */

int K;                            /* number of sets to partition numbers into */
int N;                      /* problem size: number of values in number array */
int TRIALS;                                        /* number of trials to run */
long long seed;                                 /* current random number seed */

long long sorted[MAXN];        /* original numbers sorted in decreasing order */
long long sumall;                              /* sum of all original numbers */
long long diffthresh;       /* difference threshold at which to terminate CKK */

long long calls[MAXK+1];   /* counters for number of partitionings for each K */

typedef struct element {           /* structure consists of a set and its sum */
  long long sum;
  long long set;
} SetArray;

SetArray tsets[MAXSETS];         /* temporary array of sets for merge sorting */
long long tsums[MAXSETS];        /* temporary array of sums for merge sorting */

int nextsum;              /* index of next empty position in subset sum array */

typedef struct heapelement{  /* heap of sums from a and b in increasing order */
  int x;                                         /* index of element in asums */
  int y;                                         /* index of element in bsums */
  long long sum;}                                      /* sum of the elements */
 Heap;

/* global data structures used by TWOPART */

long long asums[MAXSETS];              /* subset sums from numbers in A array */
long long bsums[MAXSETS];              /* subset sums from numbers in B array */
long long csums[MAXSETS];              /* subset sums from numbers in C array */
long long dsums[MAXSETS];              /* subset sums from numbers in D array */
int nasums, nbsums, ncsums, ndsums; /* number of sums in each of above arrays */
Heap abheap[MAXSETS], cdheap[MAXSETS]; /*heap of subset sums from pairs of sums*/
int abheapsize, cdheapsize;                       /* size of respective heaps */
long long alpha;          /* best difference found for 2-way CKK partitioning */

/* global data structures used by CGA */

int crossover[MAXK+1] =   /* largest N for which CGA or CKK is faster than SS */
  {0,0,11,10,12,14,15,17,18,20,22};

long long A[MAXN]; /* global array of numbers to be partitioned by CKK or CGA */
long long sums[MAXK];                                 /* array of subset sums */
int cgan;                              /* number of numbers to be partitioned */
int cgak;                         /* number of subsets to partition them into */
long long beta;                                 /* best solution found by CGA */
long long alpha; /* upper bound, which is maximum subset sum of completed subsets*/

/* INIT initializes it's argument array with random numbers from zero to
   2^{MAXNUM}-1.  It returns their sum. */

long long init (long long a[MAXN], int n)

{int i;                                                   /* index into array */
 long long sum;                                         /* sum of all numbers */

 sum = 0;                                    /* initialize sum of all numbers */
 for (i = 0; i < n; i++)                         /* for each element of array */
   {seed = (ACONST * seed + C) & MASK;        /* next seed in random sequence */
     a[i] = seed;                       /* random value from zero to 2^{48}-1 */

     sum += a[i];}                              /* compute sum of all numbers */
  return (sum);}                                      /* return sum of numbers */

 /* SORTDOWN takes an array of longs, and the length of the array,
    and sorts the array in decreasing order using insertion sort. */

 void sortdown (long long a[MAXN], int n)

 {int i,j;                                   /* indices into array for sorting */
  long long temp;                              /* temporary value for swapping */

  for (i = 1; i < n; i++)
    {temp = a[i];
      for (j = i - 1; j >= 0; j--)
        if (a[j] < temp) a[j+1] = a[j];
        else break;
      a[j+1] = temp;}}

 /* SORTUP takes an array of longs, and the length of the array,
    and sorts the array in increasing order using insertion sort. */

 void sortsup (long long a[MAXN], int n)

 {int i,j;                                   /* indices into array for sorting */
  long long temp;                              /* temporary value for swapping */

  for (i = 1; i < n; i++)
    {temp = a[i];
      for (j = i - 1; j >= 0; j--)
        if (a[j] > temp) a[j+1] = a[j];
        else break;
      a[j+1] = temp;}}

 /* INSERT takes an array index FIRST, a new subset sum VECTOR, and an array A of
    vectors sorted in decreasing order of largest elements, and modifies array A
    by inserting the new vector in sorted order by largest element. */

 void insert (int first, int n, int k, long long vector[MAXK], long long vects[MAXN][MAXK])

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

 long long kk (long long a[MAXN], int n, int k, long long sum)

 {long long vects[MAXN][MAXK];           /* array of vectors N long and K wide */
  long long vector[MAXK];                  /* combination of first two vectors */
  int i, j;                                                 /* utility indices */
  long long max;                       /* absolute value of maximum subset sum */

  if (k >= n) return ((long long) a[n-1]); /*largest subset sum is largest number*/

  for (i = 0; i < n; i++)                   /* create initial array of vectors */
    {vects[i][0] = a[i];                                    /* copy to vectors */
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

 /* KK2 takes an array of numbers sorted in decreasing order, and finds the
    two-way KK partition from that array, returning its subset sum difference. It
    iteratively takes the first two numbers in the list, and replaces them with
    their difference, inserted in sorted order in the remaining list. */

long long kk2 (a, n)

     long long a[MAXN];                        /* array of subset sum vectors */
     int n;                                    /* number of elements in array */

{long long b[MAXN];                      /* copy of original array to destroy */
  int i;                                                  /* index to vectors */
  int j;                         /* index to individual subset sums in vector */
  long long diffOneTwo;  /* new number which is difference between first two numbers */

  for (i = 0; i < n; i++)      /* copy original array so as not to destroy it */
    b[i] = a[i];
 for (i = 0; i < n-1; i++)                 /* combine largest numbers in turn */
   {diffOneTwo = b[i] - b[i+1];           /* new number is difference of two largest */
   for (j = i+2; j < n ; j++)                    /* for each remaining number */
     if (diffOneTwo >= b[j]) break;            /* found correct place for new number */
     else b[j-1] = b[j];   /* new number is smaller than next, move next back */
   b[j-1] = diffOneTwo;}                     /* place new number in correct position */
 return (b[n-1]);}      /* one number left, its value is partition difference */

/* GENSETS takes an array of integers, an array of subsets to fill, the first
   index and last indices into the integer array, the sum of the numbers
   included so far, and the characteristic function of the elements included so
   far. It generates all subsets of the remaining numbers in the NUMbers array
   from index FIRST to index LAST, storing them in the SET array.  It places the
   numbers in consecutive locations indexed by the global variable NEXTSUM, leaving
   in NEXTSUM the number of sets generated. */

void gensets (long long num[MAXN],
              SetArray set[MAXSETS],
              int first, int last,
              long long cursum, long long curset)

{if (first > last)                                        /* set is completed */
  {set[nextsum].sum = cursum;     /* store subset sum in array of subset sums */
  set[nextsum++].set = curset;}       /* store characteristic function of set */
 else                                                /* set not yet completed */
   {gensets (num, set, first+1, last, cursum, curset);     /* exclude next element */
   gensets (num, set, first+1, last, cursum + num[first], curset | ((long long) 1 << first));}}

/* SORTSETS sorts an array of subset sums, and the corresponding subsets
   themselves, in increasing order using mergesort.  It takes the indices of the
   FIRST and LAST elements to be sorted, and sorts those elements in increasing
   order.  It uses a temporary array to do the merge step.  */

void sortsets (SetArray sets[MAXSETS], int first, int last)

{int middle;                                      /* index to middle of array */
 int head1, head2; /* indices of current heads of sorted lists for merge step */
 int tail;                       /* index of next available position in tempa */
 SetArray temp;                 /* temporary pair of set amd sum for swapping */
 int i, j;                              /* utility indices for insertion sort */

 if (last - first <= 50)                      /* small number of sets to sort */
   {for (i = first+1; i <= last; i++)             /* from 2nd set to last set */
       {temp = sets[i];                            /* save set to be inserted */
         for (j = i - 1; j >= first; j--) /* find position for set in sorted order */
           if (sets[j].sum > temp.sum) sets[j+1] = sets[j]; /*move larger sets forward*/
           else break;                 /* found correct position for next set */
         sets[j+1] = temp;}                       /* insert into sorted order */
     return;}

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

void sortsums (long long sums[MAXSETS], int first, int last)

{int middle;                                      /* index to middle of array */
 int head1, head2; /* indices of current heads of sorted lists for merge step */
 int tail;                       /* index of next available position in tempa */
 long long temp;                /* temporary pair of set amd sum for swapping */
 int i, j;                              /* utility indices for insertion sort */

 if (last - first <= 50)                      /* small number of sums to sort */
   {for (i = first+1; i <= last; i++)             /* from 2nd sum to last sum */
       {temp = sums[i];                            /* save sum to be inserted */
         for (j = i - 1; j >= first; j--) /* find position for set in sorted order */
           if (sums[j] > temp) sums[j+1] = sums[j]; /*move larger sets forward*/
           else break;                 /* found correct position for next set */
         sums[j+1] = temp;}                       /* insert into sorted order */
     return;}

 middle = (first + last) / 2;             /* compute middle point of subarray */
 sortsums (sums, first, middle);         /* recursively merge sort first half */
 sortsums (sums, middle+1, last);       /* recursively merge sort second half */

 head1 = first;                                         /* head of first list */
 head2 = middle+1;                                     /* head of second list */
 tail = first;               /* first empty position in sorted temporary list */

 while (head1 <= middle && head2 <= last)    /* neither of the lists is empty */
   if (sums[head1] <= sums[head2])           /* head of first list is smaller */
     tsums[tail++] = sums[head1++]; /* grab smaller element and advance head1 */
   else
     tsums[tail++] = sums[head2++]; /* grab smaller element and advance head2 */

 if (head2 > last)                         /* second list was exhausted first */
   while (head1 <= middle)                   /* until first list is exhausted */
     tsums[tail++] = sums[head1++];        /* copy first list to end of merge */

 head1 = first;                        /* reset pointer to head of both lists */
 while (head1 < head2)                            /* until end of second list */
   {sums[head1] = tsums[head1]; head1++;}} /* copy merged list back into original array*/

/* GENSUMS takes an array X of numbers, the index of the LAST number, and
   generates an array SUMS of all sums that can be achieved by adding together
   numbers from X.  Each number in X can only be used at most once.  To allow
   the recursion, additional arguments include an index NEXTX to the next
   element of the array, and the SUMSOFAR achieved down this path.  There is
   also a global pointer NEXTSUM to the next empty element of the SUMS array. At
   the end, it is equal to the number of sums created. */

void gensums (long long nums[MAXN],                          /* array of original numbers */
long long sums[MAXSETS],                                   /* array of sums */
int nextx,                            /* pointer to next element of X array */
long long sumsofar,             /* sum so far of elements in current subset */
int last)                                 /* index of last element in array */

{if (nextx == last)                          /* reached last element of array */
   {sums[nextsum++] = sumsofar;                     /* don't add last element */
    sums[nextsum++] = sumsofar + nums[nextx];}        /* add last element to sum */
  else {gensums (nums, sums, nextx+1, sumsofar, last);     /* don't add last element */
    gensums (nums, sums, nextx+1, sumsofar + nums[nextx], last);}} /* add last element */

/* CKK takes an array of numbers A sorted in decreasing order, its
   length N, and their TOTAL sum, and finds the best partition,
   leaving the resulting difference in the global variable ALPHA. It
   runs branch-and-bound, starting with the Karmarkar-Karp
   solution. At each point, the largest two remaining numbers are
   selected, and replaced with either their difference or their sum,
   representing assigning them to different sets or the same set,
   respectively. */

void ckk (long long a[MAXN],                                        /* array of numbers */
int n,                                         /* number of elements in array */
long long total)                                        /* sum of all numbers */

{long long diff2;                          /* difference of largest 2 numbers */
 long long sum2;                                  /* sum of largest 2 numbers */
 long long difference;                             /* difference of partition */
 long long b[MAXN];                   /* new copy of list for recursive calls */
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
    if (alpha <= diffthresh) return;} /*good enough solution based on maxsofar*/

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

/* TWOPART takes an array of NUMberS, the Number of numbers in the
   array, their SUM, the maximum subset sum in the previously
   completed subsets (MAXSOFAR) in the current solution path, and the
   largest subset sum in the best complete solution found so far
   (BESTSOFAR), and partitions them two ways using the Schroeppel and
   Shamir algorithm.  It returns the larger subset sum in the first
   partition it finds whose larger subset sum is less than or equal to
   MAXSOFAR. */

long long twopart
(long long nums[MAXN],                 /* array of integers to be partitioned */
 int n,                               /* number of integers to be partitioned */
 long long sum,                                    /* sum of all the integers */
 long long maxsofar, /* largest subset sum so far in the current solution path*/
 long long bestsofar) /* largest subset sum in best complete solution found so far */

{int i;                                                      /* utility index */
 long long best;                           /* best partition difference found */
 long long target1;       /* two target values, equal if sum is even */
 int dpointer;                 /* index to largest number in dsums <= target2 */
 int nexty;           /* temporary for y component of next element of CD heap */
 int child, parent;              /* indices of parent and child nodes in heap */
 Heap next;                                     /* new element to add to heap */
 long long newsum;                         /* new subset sum being considered */
 long long compsum;                           /* subset sum of complement set */
 long long larger;                            /* larger of newsum and compsum */

 calls[2]++;                              /* count number of 2-way partitions */

 if (n <= crossover[2])    /* small number of elements, run CKK instead of SS */
   {alpha = MASK;                         /* initial value of best difference */
   for (i = 0; i < n; i++) /* */
     A[i] = nums[i];                        /* copy numbers into global array */
   if (n < 5) alpha = kk2 (A, n);    /* ckk only works with 5 or more numbers */
   else
     {diffthresh = 2 * maxsofar - sum;      /* larger sum will be <= maxsofar */
       ckk (A, n, sum);}   /* run CKK to find best 2-way partition difference */
   return ((sum + alpha) / 2);}  /* return larger subset of optimal partition */

 best = (sum + kk2(nums, n)) / 2;                /* larger set in KK solution */
 if (best <= maxsofar) return (best);           /* KK solution is good enough */
 else if (best > bestsofar) best = bestsofar;

 target1 = sum / 2;                                   /* smaller target value */
 nums = &nums[1];                      /* exclude largest number from subsets */
 n = n - 1;                                              /*  one fewer number */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, dsums, 0, 0ll, n-1-n*3/4);        /* generate sums of elements */
 ndsums = nextsum;                                /* number of sums generated */
 sortsums (dsums, 0, ndsums-1);            /* merge sort array of sums from d */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, csums, n-n*3/4, 0ll, n-1-n/2);    /* generate sums of elements */
 ncsums = nextsum;                                /* number of sums generated */
 sortsums (csums, 0, ncsums-1);            /* merge sort array of sums from c */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, bsums, n-n/2, 0ll, n-1-n/4);      /* generate sums of elements */
 nbsums = nextsum;                                /* number of sums generated */
 sortsums (bsums, 0, nbsums-1);            /* merge sort array of sums from b */

 nextsum = 0;                               /* initialize index to sums array */
 gensums(nums, asums, n-n/4, 0ll, n-1);          /* generate sums of elements */
 nasums = nextsum;                                /* number of sums generated */
 sortsums (asums, 0, nasums-1);            /* merge sort array of sums from a */

 for (i = 0; i < nasums; i++)                                 /* make ab heap */
   {abheap[i].x = i;                         /* aptr pointer is same as index */
    abheap[i].y = 0;              /* smallest value of bsums is first element */
    abheap[i].sum = asums[i];        /* asums[x] + bsums[0], but bsums[0] = 0 */
    if (abheap[i].sum >= best) break;}       /* sum too big, done making heap */
 abheapsize = i;                      /* number of elements in heap initially */
                /* make CD maxheap, only with elements <= larger target value */
 dpointer = ndsums - 1;                    /* index to largest value in dsums */
 cdheapsize = 0;                                   /* no elements in heap yet */
 for (i = 0; i < ncsums; i++)            /* process csums in increasing order */
   {if (csums[i] >= best) break;            /* no more entries to add to heap */
    else cdheapsize++;                /* index of next empty location in heap */
    next.x = i;                                          /* cpointer is index */
    while (csums[i] + dsums[dpointer] >= best) dpointer--;    /* sums too big */
    next.y = dpointer;                    /* largest value of dsums <= target */
    next.sum = csums[i] + dsums[dpointer];         /* sum of new heap element */
    child = cdheapsize-1;                       /* index of last node in heap */
    while (child > 0)                  /* bubble child up to correct location */
      {parent = (child - 1) / 2;          /* index in heap of parent of child */
       if (next.sum > cdheap[parent].sum)   /* new element is too low in heap */
         {cdheap[child] = cdheap[parent];                 /* move parent down */
          child = parent;}                        /* parent is new child node */
       else break;}             /* found correct location of new heap element */
    cdheap[child] = next;}                             /* insert it into heap */

 while (abheapsize > 0 && cdheapsize > 0)  /* there are still valid sums left */
   {newsum = abheap[0].sum + cdheap[0].sum; /* next ab number + next cd number*/
    compsum = sum - newsum;                   /* subset sum of complement set */
    if (compsum >= newsum) larger = compsum;    /*  complement has larger sum */
    else larger = newsum;                      /* original set has larger sum */
    if (larger < best)                               /* better solution found */
      {best = larger;                                   /* update best so far */
        if (larger <= maxsofar) return (larger);} /*good enough partition found*/
    if (newsum < target1)               /* construct next element for AB heap */
      {if (abheap[0].y + 1 < nbsums)            /* still more b pointers left */
          {next.x = abheap[0].x;                  /* x pointer stays the same */
           next.y = abheap[0].y + 1;          /* set to next larger b pointer */
           next.sum = asums[next.x] + bsums[next.y];} /*sum of new heap element*/
       else                                  /* remove this element from heap */
         {abheapsize--;                        /* heap has one fewer elements */
          next = abheap[abheapsize];} /*replace top with last element of heap */
                       /* percolate new top of heap down to its correct place */
        parent = 0;
        while (2 * parent + 1 < abheapsize)  /* still at least one child left */
          {if (2 * parent + 2 == abheapsize)       /* node only has one child */
              child = 2 * parent + 1;
           else if (abheap[2*parent+1].sum < abheap[2*parent+2].sum) child = 2*parent+1;
           else child = 2 * parent + 2;                 /* find smaller child */
           if (next.sum <= abheap[child].sum) break;           /* found place */
           else
             {abheap[parent] = abheap[child];    /* bubble child up to parent */
              parent = child;}}                    /* child is now new parent */
        abheap[parent] = next;}     /* insert next element in parent position */

    else                 /* sum is greater than target, remove top of cd heap */
      {nexty = cdheap[0].y-1;
       while (nexty >= 0 &&                      /* still sums in this column */
              csums[cdheap[0].x] + dsums[nexty] + abheap[0].sum >= best)
         nexty--;                               /* total sum is still too big */
       if (nexty >= 0)                     /* didn't run out of sums in colun */
         {next.y = nexty;             /* set y coordinate in new heap element */
          next.x = cdheap[0].x;                      /* x pointer is the same */
          next.sum = csums[next.x] + dsums[nexty];}    /* sum of heap element */
       else                                  /* remove this element from heap */
         {cdheapsize--;                         /* heap has one fewer element */
          next = cdheap[cdheapsize];}    /* make new top of heap last element */
                         /* sink new top of heap down to its correct position */
        parent = 0;
        while (2 * parent + 1 < cdheapsize)  /* still at least one child left */
          {if (2 * parent + 2 == cdheapsize)       /* node only has one child */
              child = 2 * parent + 1;
           else if (cdheap[2*parent+1].sum > cdheap[2*parent+2].sum) child = 2*parent+1;
           else child = 2 * parent + 2;                  /* find larger child */
           if (next.sum >= cdheap[child].sum) break;           /* found place */
           else
             {cdheap[parent] = cdheap[child];
              parent = child;}}                    /* child is now new parent */
        cdheap[parent] = next;}}    /* insert next element in parent position */
 return (best);}                              /* return best subset sum found */

/* CGA accesses a global vector of CGAK subset SUMS, sorted in
   increasing order, the NEXT index in the array A of numbers, sorted
   in decreasing order, and the FIRST possible position to place a new
   number, which is the last subset sum of zero in the array. It adds
   A[next] to each subset in the vector, as long as the resulting sum
   is strictly less than the best solution found so far, BETA.  For
   each new partial partition, it recursively searches all possible
   completions of the partial partition that could be better than
   BETA.  It will only place the new element is one empty subset sum,
   to avoid generating permutations of the same partition. As soon as
   it finds a partition whose largest subset sum is less than or equal
   to the global variable ALPHA, or an optimal solution, it
   returns. As a side-effect, it sets BETA to the maximum subset sum
   in the best solution found. */

void cga (int next, int first)

{long long largestOfLastTwo;                      /* largest subset sum of last two subsets */
 long long y, z;                            /* subset sums of last two subsets */
 int i,j;                                                     /* utility index */



 if (next == cgan - 2)                                /* only two numbers left */
   {z = sums[0] + A[next];   /* add next to last number to smallest subset sum */
     if (z <= sums[1]) largestOfLastTwo = z + A[next+1];  /* add last number to same subset */
     else
       {y = sums[1] + A[next+1];  /* add last number to second smallest subset */
         if (z >= y) largestOfLastTwo = z;
         else largestOfLastTwo = y;}
     if (largestOfLastTwo <= sums[cgak-1]) largestOfLastTwo = sums[cgak-1];    /* new largest subset sum */
     if (largestOfLastTwo < beta) beta = largestOfLastTwo;     /* found better solution than best so far */
     return;}                             /* this solution is complete, return */

 for (i = first; i < cgak; i++) /* consider placing next number in every subset*/
   {largestOfLastTwo = A[next] + sums[i];              /* add next number to current subset */
     if (largestOfLastTwo < beta)            /* new subset sum is smaller than smallest max */
       {for (j = i; j < cgak - 1; j++) /*find place for new subset in sorted order*/
           if (largestOfLastTwo > sums[j+1]) sums[j] = sums[j+1]; /*copy smaller elements back*/
           else break;                       /* found correct place, break out */
         sums[j] = largestOfLastTwo;              /* add new subset sum in correct position */
         if (first > 0 && i == first) cga (next + 1, first - 1); /*added to zero sum */
         else cga (next + 1, first);                  /* added to non-zero sum */
         if (beta <= alpha) return;              /* good enough solution found */
         largestOfLastTwo = sums[j] - A[next];                 /* original element of array */
         for (j = j; j > i; j--)              /* insert back into sorted order */
           sums[j] = sums[j-1];                /* move larger elements forward */
         sums[i] = largestOfLastTwo;                   /* restore original element in array */
         if (sums[cgak-1] >= beta) return;}   /* largest >= max of best so far */
     else break;} /* new element >= beta, no point in adding to larger subsets */
 return;}                          /* return, without finding optimal solution */

/* PART takes an array NUM of integers sorted in decreasing order, the
   length N of the vector, the SUM of all the numbers, a number K of
   subsets to partition it into, the largest subset sum of the
   completed subsets in the current solution path (MAXSOFAR), and the
   largest subset sum in the best complete solution found so far
   (BESTSOFAR).  It returns either the largest subset sum of the first
   solution it finds where all subset sums are less than or equal to
   MAXSOFAR, the largest subset sum of the optimal k-way partition if
   it is less than BESTSOFAR, or BESTSOFAR otherwise. */

long long
part (long long num[MAXN],              /* array of integers to be partitioned */
      int n,                                    /* number of integers in array */
      long long sum,                                /* sum of all the integers */
      int k,                       /* number of ways to partition the integers */
      long long maxsofar, /*largest sum among completed subsets in current path*/
      long long bestsofar) /* largest subset sum in best complete solution found so far */

{long long lower;                       /* lower bound of top-level subset sum */
  long long upper;                      /* upper bound of top-level subset sum */
  long long perfect;                /* maximum subset sum in perfect partition */

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

  long long x[MAXN]; /* array to contain integers excluded from current subset */
  int nx;                              /* number of integers in excluded array */
  int first;                                /* index of first subset in cdlist */
  int next;                         /* index of first empty location in cdlist */
  long long topsum;        /* sum of combination subset from top of both heaps */
  int index;                                              /* index into cdlist */
  long long sol;             /* maximum subset sum from recursive subpartition */
  long long set;           /* characteristic function of set being constructed */
  long long subsum;                              /* subset sum of first subset */
  int parent, child;              /* indices of parent and child nodes in heap */
  int i;                                                      /* utility index */
  long long best;      /* largest subset sum in best partition of argument set */
  int nexty;           /* temporary for y component of next element of CD heap */
  long long largest;                         /* largest number in original set */
  long long newmax; /* new maximum subset sum of completed subsets on current path*/
  long long sumrest;       /* sum of remaining included items in current subst */
  long long rest;                                  /* sum of excluded elements */
  long long lastex;       /* smallest excluded number so far in current subset */

  calls[k]++;                           /* count number of calls at this level */

  if (n <= k) return(num[0]); /* no more numbers than subsets, return largest */

  if (n <= crossover[k])             /* use CGA instead of recursive algorithm */
    {beta = bestsofar; /*initial upper bound is max sum in best solution so far*/
      for (i = 1; i < n; i++)            /* largest number is never referenced */
        A[i] = num[i];                       /* copy numbers into global array */
      cgan = n;         /* number of integers to be partitioned for efficiency */
      cgak = k;               /* number of subset into which to partition them */
      alpha = maxsofar; /* maximum subset sum of completed subsets in current solution */
      for (i = 0; i < k-1; i++)                       /* initialize sums array */
        sums[i] = 0;
      sums[k-1] = num[0];              /* put largest number in largest subset */
      cga (1, k-2);            /* find solution with complete greedy algorithm */
      return (beta);}          /* return largest subset in best solution found */

  best = kk (num, n, k, sum);     /* Karmarkar-Karp solution is initial approx */
  if (n <= k + 2) return (best); /*kk is optimal if k <= n + 2. ONLY NEEDED FOR DEBUGGING*/
  if (sum % k == 0) perfect = sum / k;                /* largest subset sum of */
  else perfect = sum / k + 1;                             /* perfect partition */
  if (num[0] > perfect) perfect = num[0];         /* can't be < largest number */
  if (perfect > maxsofar) maxsofar = perfect;  /* can't do better than perfect */
  if (best <= maxsofar) return (best);     /* kk solution is perfect partition */
  if (best >= bestsofar) best = bestsofar; /* don't look for solutions worse than best */
  else bestsofar = best;                         /* kk is best solution so far */

  largest = num[0];                                   /* largest number in set */
  num = &num[1];                                        /* skip largest number */
  n = n - 1;                      /* exclude largest number from consideration */
  upper = best - 1 - largest;               /* upper bound on first subset sum */
  lower = sum - (k-1) * (best - 1) - largest; /*lower bound on first subset sum*/

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets(num, d, 0, n-1-n*3/4, 0ll, 0ll);         /* generate sets of elements */
  ndsets = nextsum;                             /* number of subsets generated */
  sortsets (d, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets(num, c, n-n*3/4, n-1-n/2, 0ll, 0ll);     /* generate sets of elements */
  ncsets = nextsum;                             /* number of subsets generated */
  sortsets (c, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets(num, b, n-n/2, n-1-n/4, 0ll, 0ll);       /* generate sets of elements */
  nbsets = nextsum;                             /* number of subsets generated */
  sortsets (b, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

  nextsum = 0;                    /* starting index to place sets in set array */
  gensets(num, a, n-n/4, n-1, 0ll, 0ll);           /* generate sets of elements */
  nasets = nextsum;                             /* number of subsets generated */
  sortsets (a, 0, nextsum - 1); /* sort subsets in increasing subset sum order */

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
      newElement.x = index;                                     /* C pointer is index */
      while (c[index].sum + d[dpointer].sum > upper) dpointer--; /* find sums <= upper */
      newElement.y = dpointer;                    /* largest value of dsums <= target */
      newElement.sum = c[index].sum + d[dpointer].sum;         /* sum of heap element */
      child = cdheapsize-1;                     /* index in heap of child node */
      while (child > 0)                 /* bubble child up to correct location */
        {parent = (child - 1) / 2;         /* index in heap of parent of child */
          if (newElement.sum > cdheap[parent].sum)        /* child is too low in heap */
            {cdheap[child] = cdheap[parent];
              child = parent;}                     /* parent is new child node */
          else break;}  /* child is in correct place, add next element to heap */
      cdheap[child] = newElement;}                             /* insert it into heap */

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
      while (cdheapsize > 0 && topsum >= lower)       /* new set within bounds */
        {if (topsum <= upper)                         /* new set within bounds */
            {cdlist[next] = cdheap[0];             /* add heap element to list */
              next = (next + 1) % MAXLIST;     /* add new set to end of ablist */
              cdlistsize++;
              if (cdlistsize > MAXLIST)
                {printf ("error: subset list size exceeded\n"); exit(0);}}

          nexty = cdheap[0].y-1;                     /* next element in column */
          while (nexty >= 0 &&                    /* still sums in this column */
              c[cdheap[0].x].sum + d[nexty].sum + abheap[0].sum > upper)
            nexty--;                             /* total sum is still too big */
          if (nexty >= 0)                   /* didn't run out of sums in colun */
            {newElement.y = nexty;            /* set y coordinate in new heap element */
             newElement.x = cdheap[0].x;                     /* x pointer is the same */
             newElement.sum = c[newElement.x].sum + d[nexty].sum;}    /* sum of heap element */
          else                              /* no more elements in this column */
            {cdheapsize--;                      /* size of heap shrinks by one */
              if (cdheapsize > 0)                   /* CDheap is not yet empty */
                newElement = cdheap[cdheapsize];/*replace top of heap with last member*/
              else break;}         /* CDheap is empty, break out of while loop */
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
                 cdheap[parent] = newElement; /*insert next element in parent position*/
          topsum = cdheap[0].sum + abheap[0].sum;}  /* sum of two top elements of heaps */

      if (cdlistsize > 0)/*any subsets left in the list between lower and upper*/
        {index = first;                                /* first subset in list */
          while (index != next)        /* there are still more subsets on list */
            {subsum = cdlist[index].sum + abheap[0].sum;   /* first subset sum */
              if (subsum >= lower && subsum <= upper)/*subset sum within bounds*/
                {set = c[cdlist[index].x].set |
                       d[cdlist[index].y].set |
                       a[abheap[0].x].set |
                       b[abheap[0].y].set;
                  nx = 0;      /* constuct array of numbers excluded from set */
                  sumrest = subsum; /* sum of all included integers except largest*/
                  lastex = sum; /* last excluded number, initially large value*/
                  for (i = 0; i < n; i++)
                    {if ((set & 1) == 0)     /* number is not in first subset */
                        {if ((num[i] >= sumrest) &&
                             (largest + subsum + num[i] - sumrest <= maxsofar)) break;
                          x[nx] = num[i];       /* add number to excluded set */
                          nx++;              /* count numbers in original set */
                          lastex = num[i];} /* update last excluded number */
                      else                    /* number is included in subset */
                        {if (largest + subsum - num[i] + lastex <= maxsofar) break;
                          sumrest -= num[i];} /* sum of remaining included elements */
                      set = set >> 1;}   /* shift to get bit for next element */
                  if (i < n)                       /* subset should be pruned */
                    {index = (index + 1) % MAXLIST;  /* next subset on CDlist */
                      continue;}                          /* skip this subset */

                  if (x[0] < bestsofar && subsum + largest + x[nx-1] > maxsofar)
                    {rest = sum-largest-subsum;  /* sum of remaining elements */
                      if (subsum + largest >= maxsofar) newmax = subsum + largest;
                      else newmax = maxsofar; /* new largest subset sum so far*/
                      if (rest / (k-1) > newmax) newmax = rest / (k-1); /*perfect part*/
                      if (k-1 == 2) sol = twopart (x, nx, rest, newmax, bestsofar);
                      else sol = part (x, nx, rest, k-1, newmax, bestsofar);
                      if (sol < bestsofar)          /* found better partition */
                        {if (sol < subsum + largest) sol = subsum + largest; /* larger sum*/
                          best = sol;     /* update best local solution */
                          if (best <= maxsofar) return (best);
                          upper = best - 1 - largest;  /* new upper bound */
                          lower = sum - (k-1) * (best - 1) - largest; /* new lb */
                          if (lower > upper) return (best); /* can't find better solution*/
                          else bestsofar = best;}}} /* update best global solution */

              index = (index + 1) % MAXLIST;}}                /* next subset on CDlist */

      if (abheap[0].y < nbsets - 1)         /* there is another element in this column */
        {newElement.x = abheap[0].x;          /* replace with next combination in this column */
          newElement.y = abheap[0].y + 1;
          newElement.sum = a[newElement.x].sum + b[newElement.y].sum;} /* sum of new subset */
      else
        {abheapsize--;                           /* decrease size of heap by one */
          if (abheapsize > 0)                           /* ABheap is not yet empty */
            newElement = abheap[abheapsize];     /* replace top of heap with last member */
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

/* This is the main function.  It calls init to generate the random numbers,
   then sorts them in decreasing order.  Next it calls KK5 to generate an
   approximate 5-way partition.  If there are more that 7 numbers, or KK5
   doesn't return an optimal partition, then it calls five-way to optimally
   partition the numbers. */

int main (int argc, char *argv[])

{long long totalmax;              /* total of all 5-way partition differences */
  long long totalcalls[MAXK+1]; /* total number of recursive calls for each k */
  int i;                                                     /* utility index */
  int trial;                                       /* number of current trial */
  long long perfect;               /* largest subset sum in perfect partition */
  long long minmax;          /* the largest subset sum in an optimal solution */

 sscanf(argv[1], "%d", &K);       /* read number of subsets from command line */
 sscanf(argv[2], "%d", &N);        /* read number of values from command line */
 sscanf(argv[3], "%d", &TRIALS);   /* read number of values from command line */

 seed = INITSEED;               /* initialize seed of random number generator */
 for (i = 2; i <= K; i++)                  /* initialize total calls counters */
   totalcalls[i] = 0;
 totalmax = 0;                      /* initialize total solution cost counter */

 for (trial = 1; trial <= TRIALS; trial++)      /* for each independent trial */
   {sumall = init (sorted, N);    /* generate random values, return their sum */
    sortdown (sorted, N);                 /* sort numbers in decreasing order */

   for (i = 2; i <= K; i++)               /* initialize total calls counters */
     calls[i] = 0;

   perfect = (sumall+(K-1)) / K;
   if (K == 2) minmax = twopart (sorted, N, sumall, sumall - sumall/2, sumall); /* optimal 2-way */
   else minmax = part (sorted, N, sumall, K, perfect, sumall); /* optimally partition K ways*/

   totalmax += minmax;              /* increment total solution cost counter */
   for (i = 2; i <= K; i++)                /* increment total calls counters */
     totalcalls[i] += calls[i];

   printf ("%3d %lld\n", trial, minmax);} /* output each individual instance */
 printf("%d %2d %lld ", K, N, totalmax);       /* solution statistics for run */
 for (i = K; i >= 2; i--)                 /* print number of calls for each k */
   printf("%lld ", totalcalls[i]);
 printf ("\n");
return 0;
}
