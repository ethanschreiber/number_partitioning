#include <stdlib.h>
#include "../../PackingUtils.hpp"
#include "IECompletionGenerator.hpp"
#include <stdint.h>
#include <iostream>
#include <vector>
#include <deque>
#include <stack>

using namespace std;
int main() {
  const int NUM_ELEMENTS=5;
  const int64_t elements[NUM_ELEMENTS] = {10,9,8,7,6};
  const int64_t LOWER_BOUND = 6;
  const int64_t UPPER_BOUND = 22;


  int64_t sum = 0;
  for (int i=0;i<NUM_ELEMENTS;i++) {
    sum += elements[i];
  }
  ss::IECompletionGenerator generator (LOWER_BOUND,UPPER_BOUND,elements,NUM_ELEMENTS,sum);

  vector<SetNodeVector> sets2;
  SetNodeVector node;

  while (generator.next(node)) {
    sets2.push_back(node);
  }
  std::sort(sets2.begin(),sets2.end(),SetNodeCardinalityComparator());

  for (size_t i=0;i<sets2.size();i++) {
    cout << "--> " << sets2[i].toString() << endl;
  }

}
