#include <numeric> // For accumulate
#include <deque>
#include <stdint.h>
#include <boost/dynamic_bitset.hpp>
#include "CachedSetsBuffered.hpp"
#include "PartitionUtils.hpp"

using ss::SetNodeBitset;
using partition::SetVector;
using partition::SetVectorIt;
using partition::SetDeque;
using partition::SetIteratorPair;

class SubsetSelector {
private:

  uint64_t m_ub;
  uint64_t m_K;
  const SetVector &m_subsets;
  vector<ss::DynamicBitset> m_isDisjoint;
  size_t m_count;
  const int LOW_CARD;
  const int HIGH_CARD;
public:
  size_t count() const {
    return m_count;
  }
  SubsetSelector(uint64_t ub, uint64_t K, const SetVector &subsets, const int lowCard)
      : m_ub(ub), m_K(K), m_subsets(subsets), m_count(0), LOW_CARD(lowCard), HIGH_CARD(lowCard + 1) {
    SimpleTimer timer;
    // Setup Disjoint Vector
    for (size_t i = 0; i < subsets.size(); i++) { // One entry for each subset
      m_isDisjoint.push_back(ss::DynamicBitset(subsets.size(), 0));
    }
    for (size_t i = 0; i < subsets.size(); i++) { // One entry for each subset
      for (size_t j = i + 1; j < subsets.size(); j++) {
        if ((subsets[i].set() & subsets[j].set()).none()) {   // If disjoint
          m_isDisjoint[i].set(j, true);                        // subset i is disjoint of subset j
          m_isDisjoint[j].set(i, true);                        // subset j is disjoint of subset i
        }
      }
    }

    cout << "Time TO Construct: " << timer.timeElapsed() << endl;
  }

  string subsetsToString(const SetVector &subsets,const int MAX_CARD) {
    std::ostringstream out;
    size_t counts[MAX_CARD+1];

    out << "[";
    memset(counts,0,sizeof(size_t) * (MAX_CARD+1));
    for (size_t i=0;i<subsets.size();i++) {
//      cout << (int) subsets[i].cardinality() << " ";
      counts[subsets[i].cardinality()]++;
    }

    for (int i=0;i<=MAX_CARD;i++) {
      if (counts[i]) {
        out << i << ": " << counts[i] << "  " ;
      }
    }
//    cout << endl << out.str() << endl;
    out << "]  ";
    return out.str();

  }
  void selectSubsets(const int k, const uint64_t sum, const SetNodeBitset &subset, const SetVector &subsets) {

    // Debug Start
    int nLeft = subset.set().size() - (int) subset.cardinality();
    int kLeft = m_K - k;
//    if (nLeft >= (kLeft * HIGH_CARD)) {
    if (k <= 10) {
      m_count++;
      cout << "|" << partition::spaces(k) << k << " " << (int) subset.cardinality() << " SS Left: " << subsetsToString(subsets,LOW_CARD) << "(" << kLeft << ", " << nLeft << ") "
           << (kLeft * HIGH_CARD) - nLeft << endl;
      if ( subset.cardinality() == 6) {
        exit(0);
      }
    }
    // Debug End

    for (size_t i = 0; i < subsets.size(); i++) {                         // For each subset
      uint64_t sumI = sum - subsets[i].sum();                         		// Compute new sum
      uint64_t lbI = partition::computeCMin(sumI, m_K - (k + 1), m_ub);  	// Minimum value for remaining elements

      if (lbI < m_ub) {
        SetNodeBitset newSubset(subset, subsets[i]);                      // Combine them

        SetVector newSubsets;                                             // Include all subsets which are still disjoint
        newSubsets.reserve(subsets.size());

        for (size_t j = i + 1; j < subsets.size(); j++) {

          uint64_t sumIJ = sumI - subsets[j].sum();                                     // Only add if not less than lb
          uint64_t lbIJ = partition::computeCMin(sumIJ, m_K - (k + 2), m_ub);           // Minimum value for remaining elements

          if (lbIJ < m_ub) {                                  // If lb would still be below ub AND
            if ((newSubset.set() & subsets[j].set()).none()) {  // If disjoint
              newSubsets.push_back(subsets[j]);                 // Add
            }
          }
        }
        selectSubsets(k + 1, sumI, newSubset, newSubsets);
      }
    }
  }

//  void selectSubsets(const int depth, const uint64_t sum, const SetNodeBitset &subset, const ss::DynamicBitset &disjointSubsets) {
//    // Debug Start
//    int nLeft = subset.set().size() - (int) subset.cardinality();
//    int kLeft = m_K - depth;
////    if ((nLeft / kLeft) > 6) {
//    if (depth <= 2) {
//      m_count++;
//      cout << "|" << partition::spaces(depth) << depth << " " << (int) subset.cardinality() << " SS Left: " << disjointSubsets.count() << "(" << kLeft << ", "
//           << nLeft << ") " << endl;
//    }
////    }
//    // Debug End
//
//    size_t i = disjointSubsets.find_first();                   // First set
//
//    while (i != ss::DynamicBitset::npos) {       // While there are disjoint subsets left
//      uint64_t sumI = sum - m_subsets[i].sum();             // Compute new sum
//      uint64_t lbI = partition::computeCMin(sumI, m_K - (depth + 1), m_ub);             // Minimum value for remaining elements
//
//      if (lbI < m_ub) {
//        const SetNodeBitset newSubset(subset, m_subsets[i]);           // Combine them
//        const ss::DynamicBitset newDisjointSubsets = disjointSubsets & m_isDisjoint[i];
//
//        selectSubsets(depth + 1, sumI, newSubset, newDisjointSubsets);
//      }
//
//      i = disjointSubsets.find_next(i);
//    }
//  }
};

int main(int argc, char *argv[]) {
  const int N = 60;
  const int K = 10;
  const int LOW_CARD = 6;
  const int NUM_SETS = 3000;
  const uint64_t S[N] = { 269371976871218, 260786717900045, 242374462877153, 240328834836393, 237512738170429, 229437967113590, 224826668205860,
                          222386574801425, 213814114975560, 213350595726533, 207667799258467, 206212442942956, 193466919245191, 186878042962176,
                          183163933391052, 177956397673062, 171243073288138, 170210829186182, 162447845321538, 159165202205751, 157952248550395,
                          151061763370712, 150490677936643, 150124824069776, 140016007861764, 139702721835791, 138217487076690, 135325785521278,
                          129097006853268, 118469058978227, 114932361670064, 110946638862368, 110188590872511, 109679932364371, 104043563803240,
                          100105851762469, 99473552299415, 98294412649532, 91374850364258, 90382443677626, 86930750925397, 77618436305273, 76576498681259,
                          69002486940939, 68634971600750, 64101275900758, 63459886835124, 62138814596074, 56940174459161, 55965142520156, 48269979061210,
                          41966776361567, 41601943234221, 38808442460149, 38685990807517, 31950913910193, 17598934570590, 9939222116071, 9154903641929,
                          7800734460763 };

  uint64_t sum = 0;
  sum = std::accumulate(S, S + N, sum);
  uint64_t upperBound = kk(S, N, K, sum);

  SimpleTimer timer;
  cout << "Building Cache...";
  cout.flush();
  partition::CachedSetsBuffered *buf = new partition::CachedSetsBufferedLowCardinalityHS(S, N, K, upperBound, NUM_SETS, LOW_CARD);

  cout << timer.timeElapsed() << " seconds." << endl;
  size_t counts[LOW_CARD + 1];
  memset(counts, 0, sizeof(size_t) * (LOW_CARD + 1));
  SetDeque subsetsDeque[LOW_CARD + 1];
//	SetDeque subsetsDeque;

  for (size_t i = 0; i < NUM_SETS; i++) {

    if (buf->next()) {
      int cardinality = buf->getMax()->cardinality();
      counts[cardinality]++;
      subsetsDeque[cardinality].push_front(*(buf->getMax()));
//		  subsetsDeque.push_front(*(buf->getMax()));

      SetIteratorPair minRange = buf->getMinRange();
      for (SetVectorIt it = minRange.first; it != minRange.second; it++) {
        counts[it->set().count()]++;
        int cardinality = it->cardinality();
        subsetsDeque[cardinality].push_back(*it);
//	      subsetsDeque.push_back(*it);
      }
    }
    if (counts[6] >= 561) {

      break;
    }
  }
  upperBound = buf->getMax()->sum();

  for (int i = 0; i <= LOW_CARD; i++) {
    if (counts[i]) {
      cout << std::setw(2) << i << ": " << counts[i] << endl;
    }
  }

//  // Put the sets into the setVector in order of cardinality descending
//  SetVector subsetsVector;
//  for (int card=LOW_CARD;card >= 0;card--) {
//    std::copy(subsetsDeque[card].begin(), subsetsDeque[card].end(), std::back_inserter(subsetsVector));
//  }

  // Put the sets into the setVector in order of cardinality ascending
  // for each card, put in order of sum descending
  SetVector subsetsVector;
  for (int card = 0; card <= LOW_CARD; card++) {
    std::copy(subsetsDeque[card].begin(), subsetsDeque[card].end(), std::back_inserter(subsetsVector));
  }

//	std::copy(subsetsDeque.begin(), subsetsDeque.end(), std::back_inserter(subsetsVector));

  SetNodeBitset subset(0, ss::DynamicBitset(N));

  cout << "Num Sets: " << subsetsVector.size() << endl;

  const int MIN_HIGH_CARD = LOW_CARD + 1;
  for (int i = 0; i <= K; i++) {
    cout << std::setw(2) << i << ": " << std::setw(6) << N - ((K - i) * MIN_HIGH_CARD) << endl;
  }

  SubsetSelector selector(upperBound, K, subsetsVector, LOW_CARD);
  selector.selectSubsets(0, sum, subset, subsetsVector);

//  ss::DynamicBitset disjointSubsets(subsetsVector.size(),0);
//  disjointSubsets.set();    // Set all to 1 initially
//  selector.selectSubsets(0,sum, subset, disjointSubsets);

  cout << "Count: " << selector.count() << endl;
  cout << "Total Count: " << timer.timeElapsed() << endl;
}
// 0:    -10
// 1:     -3
// 2:      4
// 3:     11
// 4:     18
// 5:     25
// 6:     32
// 7:     39
// 8:     46
// 9:     53
//10:     60

