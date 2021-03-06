/*
 * MoffittCachedIE.cpp
 *
 *  Created on: Dec 24, 2013
 *      Author: ethan
 */

#include "CachedIETrees.hpp"
#include <numeric> // for accumulate


namespace partition {
// ============================================================================
//
// C a c h e d I E
//
// ============================================================================

CachedIETrees::CachedIETrees(const vector<uint64_t> &S, int K) :
      m_S(S.begin(),S.end()), m_elementsSum(std::accumulate(S.begin(),S.end(),(uint64_t) 0)), m_K(K) {

	m_roots = new CachedIENode*[S.size()];
  for (size_t i=0;i<S.size();i++) { 	// Set all roots to null initially
    m_roots[i] = NULL;
  }

  d_numSubsetsPerCard.resize(S.size());

  for (size_t i=0;i<S.size();i++) {
  	d_numSubsetsPerCard[i] = 0;
  }

}

CachedIETrees::~CachedIETrees() {
	clear();
	delete [] m_roots;
}

void CachedIETrees::push_back(const ss::SetNodeBitset &node) {

  const ss::DynamicBitset &set = node.set();    // Reference to bit set

  size_t cardinality = set.count();
  d_numSubsetsPerCard[cardinality]++;
  CachedIENode ** ptr = &m_roots[cardinality];  // pointer to root pointer based on cardinality
  size_t idx = set.find_first();                // Largest element is first

  while (idx != ss::DynamicBitset::npos){       // While there are elements left

    if (*ptr == NULL || idx < (*ptr)->ieIdx) {  // If NULL or less than idx

      *ptr = new CachedIENode(idx, NULL, *ptr); // Insert new node with exclusion pointing to old one or NULL
      ptr = &((*ptr)->included);                // Point to included
      idx = set.find_next(idx);                 // Next largest
    } else if (idx == (*ptr)->ieIdx) {          // If there is a match
      ptr = &((*ptr)->included);                // Point to inclusion
      idx = set.find_next(idx);                 // Next largest
    } else {                                    // if idx > ptr->ieIdx, we exclude
      ptr = &((*ptr)->excluded);                // Point to exclusion
    }
  }
}


void CachedIETrees::clear(CachedIENode *node) {
  if (node != NULL) {
    clear(node->included);
    clear(node->excluded);
    delete node;
  }
}

void CachedIETrees::clear() {
	for (size_t i=0;i<m_S.size();i++) {
		if (m_roots[i] != NULL) {
			clear(m_roots[i]);
			m_roots[i] = NULL;
		}
	}
}

const CachedIENode* CachedIETrees::getRoot(size_t cardinality) const {
  return m_roots[cardinality];
}

size_t CachedIETrees::getNumSubsets(size_t cardinality) {
  if (cardinality < m_S.size()) {
    return d_numSubsetsPerCard[cardinality];
  } else {
    return 0;
  }
}

// Return the minimal cardinality of any set
// stored in this data structure.
// returns N if no sets are stored
size_t CachedIETrees::getMinCardinality() const {
	for (size_t i=0;i<m_S.size();i++) {
		if (getRoot(i) != NULL) {
			return i;
		}
	}
	return m_S.size();
}

// Return the maximal cardinality of any set
// stored in this data structure.
// returns N if no sets are stored
size_t CachedIETrees::getMaxCardinality() const {
	for (int i=m_S.size()-1;i>=0;i--) {
		if (getRoot(i) != NULL) {
			return i;
		}
	}
	return 0;
}

bool CachedIETrees::addSetsToCache() {

  bool hasNext = getCachedSetsBuffered()->next();                         // Get next highest max and corresponding mins

  if (hasNext) {  // If there are more sets, add them to the CIE Tree.
//    size_t card = (*(getCachedSetsBuffered()->getMax())).set().size();
  //  if ( card >= 10 && card < 31) {
  //  	cout << "SHIT: " << card << endl;
  //  }
    push_back(*(getCachedSetsBuffered()->getMax()));         // Add new max element (sum >= perfect) to cache

    // Add the range of min elements (sum < perfect) to the cache
    SetIteratorPair minRange = getCachedSetsBuffered()->getMinRange();
    for (SetVectorIt it = minRange.first; it!= minRange.second;it++) {
      push_back(*it);
    }
  }
  return hasNext;
}


#if PAPER_DEBUG == 1
void CachedIETrees::toDot(std::ostringstream &out, string path, CachedIENode  *node) {

  std::ostringstream convert;

  string newPath = path + "_" + std::to_string(m_S[node->ieIdx]);
  out << "\n   \"" << newPath << "\" [label=" << m_S[node->ieIdx] << " shape=none ]" << endl;
  out << "   \"" << path << "\" -> \"" << newPath << "\"" << endl;

  if (node->included != NULL) {
    toDot(out,newPath,node->included);
  }
  if (node->excluded != NULL) {
    newPath += "e";
    out << "\n   \"" << newPath << "\" [label=" << m_S[node->ieIdx] << " shape=none]" << endl;
    out << "   \"" << path << "\" -> \"" << newPath << "\" [style=dashed]" << endl;
    toDot(out,newPath,node->excluded);
  }
}
string CachedIETrees::toDot() {
  std::ostringstream out;
  out << "digraph CardTrees { "<< endl;
  for (size_t i=0;i<m_S.size();i++) {
    if (m_roots[i] != NULL) {

    	std::ostringstream path;
    	path << "S" << i;

    	out << "\n   \"" << path.str() << "\" [label=\"Cardinality " << i << "\" shape=none ]" << endl;

    	toDot(out,path.str(),m_roots[i]);
    }
  }
  out << "}" << endl;
  return out.str();
}

string CachedIETrees::printSets(uint64_t min, uint64_t max)  {


	uint64_t perfect = (m_elementsSum + m_K-1) / m_K;

  vector<ss::SetNodeBitset> nodes;
  ss::generateSetsSS(&m_S[0],m_S.size(),min,max,nodes);

  std::sort(nodes.begin(),nodes.end());
  std::ostringstream out;

  size_t MAX=0;

  vector< vector<uint64_t> > numMatrix;
  vector<uint64_t> sumVector;
  vector<int> iterationVector;


  // -------------------------------
  // Fill in numMatrix and sumVector
  // -------------------------------
  size_t iterationIdx = 0;
  for (ss::SetNodeBitset node : nodes) {
  	numMatrix.push_back(vector<uint64_t>());
    ss::DynamicBitset bitset = node.set();
    size_t idx = bitset.find_first();

    uint64_t sum = 0;
    while (idx != ss::SetNodeBitset::npos){

      sum += m_S[idx];
      numMatrix.back().push_back(m_S[idx]);
      idx = bitset.find_next(idx);
    }
    sumVector.push_back(sum);

    if (sum < perfect) {
    	iterationVector.push_back(-1);
    } else {


    	if (sumVector.back() != sumVector[sumVector.size()-2]) {
    		iterationIdx++;
    	}
    	iterationVector.push_back(iterationIdx);
    }

    MAX = std::max(MAX,numMatrix.back().size());
  }

  // -----------------------
  // Fill in IterationMatrix
  // -----------------------
  int forwardIdx=0;		// For filling in iteration matrix
  int reverseIdx=0;

  while (forwardIdx < (int) iterationVector.size() && iterationVector[forwardIdx] == -1) {
  	reverseIdx = forwardIdx;
  	forwardIdx++;
  }


  while (forwardIdx < (int) iterationVector.size()) {
  	uint64_t CMin = computeCMin(m_elementsSum,m_K,sumVector[forwardIdx]+1);
  	while (reverseIdx >= 0 && sumVector[reverseIdx] >= CMin) {
  		iterationVector[reverseIdx] = iterationVector[forwardIdx];
  		reverseIdx--;
  	}
  	forwardIdx++;
  }

  // -------------------------------
  // Print Sets Simply First to cout
  // -------------------------------
  for (ss::SetNodeBitset node : nodes) {
  	cout << node.toString(&m_S[0]) << " ";
  	if (node.sum() >= perfect) {
  		cout << "Min = " << computeCMin(m_elementsSum,m_K,node.sum());
  	}
  	cout << endl;
  }
  cout << "# Sets: " << nodes.size() << endl;
  cout << endl;

  // ------------------
	// Print table header
	// ------------------
	out << "\\begin{tabular}{|ccl|}" << endl
			<< "  \\hline" << endl
			<< "  \\textbf{Sum} & \\textbf{Iter} & \\textbf{Sets}\\\\"   << endl
			<< "  \\hline" << endl;



	// --------------
	// print out rows
	// --------------

	bool foundPerfect = false;
	for (size_t i=0;i<sumVector.size();i++) {
		if (!foundPerfect && sumVector[i] >= perfect) {
			foundPerfect = true;
			out << "  \\hline" << endl;
		}
		out << "  ";
		if (i % 2 == 0) {
			out << "\\rowcolor{Gray} ";
		} else {
			out << "                ";
		}
		out << sumVector[i] << " & " << iterationVector[i] << " & ";

		vector<uint64_t> &elements = numMatrix[i];
		out << "\\{" << elements[0];

		for (size_t j=1;j<elements.size();j++) {
			out << ", " << elements[j];
		}
		out << "\\}\\\\" << endl;
	}


	// ------------------
	// Print table footer
	// ------------------
	  out << "  \\hline" << endl
	  		<< "\\end{tabular}" << endl;


//  // ------------------
//  // Print table header
//  // ------------------
//  out << "\\begin{tabular}{|r|";
//  for (size_t col=0;col<numMatrix.size();col++) {
//  	if (col % 2 == 0)
//  		out << "g";
//  	else {
//  		out << "w";
//  	}
//  }
//  out << "|}" << endl
//  		<< "  \\hline" << endl;
//
//  // ----------------------------
//  // print out sums col at a time
//  // ----------------------------
//  cout << "Sums: " << sumVector.size() << endl;
//  out << "  \\rule{0pt}{2.5ex}\\textbf{Sum}";
//
//  for (size_t col=0;col < sumVector.size();col++) {
//  	if (col == sumVector.size() -1) {
//  		out << " & " << std::setw(3) << sumVector[col];
//  	} else {
//  		out << " & \\mc{" << std::setw(3) << sumVector[col] << "}";
//  	}
//  }
//  out << "\\\\" <<endl
//  		<< "  [.2ex]\\hline" << endl;
//
//  // -----------------------------
//  // print out iters col at a time
//  // -----------------------------
//
//  out << "  \\rule{0pt}{2.5ex}\\textbf{Iter}";
//  cout << "Iter: " << iterationVector.size() << endl;
//  for (size_t col=0;col < iterationVector.size();col++) {
//
//  	if (col == iterationVector.size() -1) {
//  		out << " & " << std::setw(3) << iterationVector[col];
//  	} else {
//  		out << " & \\mc{" << std::setw(3) << iterationVector[col] << "}";
//  	}
//  }
//  out << "\\\\" <<endl
//  		<< "  [.2ex]\\hline" << endl;
//
//
//  // ----------------------------
//  // Print out sets col at a time
//  // ----------------------------
//  for (size_t row=0;row<MAX;row++) {
//
//  	if (row == 0) {
//			out << "  \\rule{0pt}{2.5ex}\\textbf{Sets}";
//		} else {
//			out << "  \\rule{0pt}{2ex}               ";
//		}
//
//  	for (size_t col = 0;col < numMatrix.size(); col++) {
//
//  		if (row < numMatrix[col].size()) {
//  			out << " & " << std::setw(3) << numMatrix[col][row];
//  		} else {
//  			out << " & " << std::setw(3) << "";
//  		}
//  	}
//  	out << "\\\\" << endl;
//  }
//  out << "  \\hline" << endl
//  		<< "\\end{tabular}" << endl;

  return out.str();
}
#endif


// ----------------------------------------------------------------------------
//
// C a c h e d I E T r e e s A l l C a r d i n a l i t y S S
//
// ----------------------------------------------------------------------------

CachedIETreesAllCardinalitySS::CachedIETreesAllCardinalitySS(const vector<uint64_t>& S, int K, uint64_t upperBound, size_t numSets)
: CachedIETrees(S,K), m_cs(new CachedSetsBufferedSS(&m_S[0], m_S.size(), m_K, upperBound, numSets)) {
}

CachedIETreesAllCardinalitySS::~CachedIETreesAllCardinalitySS() {
  delete m_cs;
}

// ----------------------------------------------------------------------------
//
// C a c h e d I E T r e e s L o w C a r d i n a l i t y H S
//
// ----------------------------------------------------------------------------
CachedIETreesLowCardinalityHS::CachedIETreesLowCardinalityHS(const vector<uint64_t>& S, const int K,
                                        const uint64_t upperBound, const size_t numSets, const size_t maxCardinality)
: CachedIETrees(S,K), m_cs(new CachedSetsBufferedLowCardinalityHS(&m_S[0], m_S.size(), m_K,upperBound,numSets,maxCardinality)){
}

CachedIETreesLowCardinalityHS::~CachedIETreesLowCardinalityHS() {
	delete m_cs;
}


} /* namespace partition */




