/*
 * vertexcover.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: fgrossmann
 */

#define OUT 1
#define DEBUG_OUT 0
#define DEBUG_OUT_2 0
#define DEBUG_OUT_3 0

#define DOM_INSERT 1
#define DOM_OUT 1
#define DOM_REMOVE 0

#define PEEK_BAG 0
#define MIN_PEEK 2000

#define TIMED 0

#include "vertex_cover_lazy.h"
#include <stdio.h>
#include <unistd.h>
#include <limits>
#include <boost/foreach.hpp>

#include <gmpxx.h>

#if TIMED
#include <time.h>
#endif

namespace sequoia {

void VertexCoverLazy::init(int argc, char **argv) {

	_myone = 1;

	int ch;
	while ((ch = getopt(argc, argv, "g:t:l:")) != -1) {
		switch (ch) {
		case 'g':
			graphFileName(optarg);
			break;
		case 't':
			treedecompositionFileName(optarg);
			break;
		default:
			std::cerr << "Error: Unknown argument: " << argv[optind]
					<< std::endl;
			break;
		}
	}

	if (argc - optind > 0) {
		std::cerr << "ERROR:  Not enough options specified" << std::endl;
	}

	std::cout << "Loading Graph from: " << graphFileName() << std::endl;

	load_graph(graphFileName());

//	std::cout << "Binary form of 2 is :" << binary(2) << std::endl;

	if (treedecompositionFileName() != NULL)
		load_treedecomposition(treedecompositionFileName());
	else
		generate_treedecomposition();

        DynProgSolver::pre_solve();
	_tab.reserve(treedecomposition()->num_vertices());
	_tabSize = new long[treedecomposition()->num_vertices()];
	_root = false;

#if PEEK_BAG
	peeked = false;
	bagSize = 0;
	maxDiff = 0;
#endif

}

void VertexCoverLazy::do_leaf(const TreeDecomposition::vertex_descriptor& t) {

//	std::cout << "LEAFING: " << t << std::endl;
	aTable *tab = new aTable();
	tabElement tE = 0;
	tabIterator *tabE3 = new tabIterator(std::make_pair(tE, 0));
	tab->insert(*tabE3);
	_tab[t] = tab;
	_tabSize[t] = 1;

}

void VertexCoverLazy::do_root(const TreeDecomposition::vertex_descriptor& t) {

//	std::cout << "ROOT: " << t << std::endl;

	std::vector<int>::iterator position1, position2;
	const Bag *tee;
	tee = treedecomposition()->bag(t);

	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!
	TreeDecomposition* endTdc = new TreeDecomposition(); //tBag.size());
	TreeDecomposition::vertex_descriptor oldVertex, nuVertex;
	std::vector<TreeDecomposition::vertex_descriptor> root_vertices;
	root_vertices.reserve(tBag.size());

	long solution = std::numeric_limits<int>::max();

	if (tBag.size() == 1) {
		std::cout << "only one element!" << std::endl;
		BOOST_FOREACH(tabIterator i, *_tab[t]) {
			std::cout << "solution array: " << i.first << " value: " << i.second
					<< std::endl;
			if (i.second < solution) {
				if (i.first == 0)
					solution = i.second;
				else
					solution = i.second + 1;
			}
		}
	} else {

		position1 = tBag.begin();
		nuVertex = endTdc->add_vertex();
		for (position2 = position1; position2 < tBag.end(); position2++) {
			endTdc->bag(nuVertex)->add(*position2);
		}
		root_vertices.push_back(nuVertex);
//		std::cout << "INSERTED initial nuVertex: "<< nuVertex << std::endl;

		_tab[nuVertex] = new aTable(*_tab[t]); // we _should_ use a new _tab, but we ain't cause we feelin dirty
		delete _tab[t];
		oldVertex = nuVertex;

//		std::cout << "ROOT BREAKPOINT 1:  " << std::endl;

		for (position1++; position1 < tBag.end(); position1++) {

			nuVertex = endTdc->add_vertex();
			for (position2 = position1; position2 < tBag.end(); position2++) {
				endTdc->bag(nuVertex)->add(*position2);
			}
			endTdc->add_edge(nuVertex, oldVertex);
			root_vertices.push_back(nuVertex);
//			std::cout << "INSERTED nuVertex: "<< nuVertex << std::endl;
			oldVertex = nuVertex;
		}

		/*
		 * NOW treedecomposition() will return the newly made root tdc
		 */
		_root_tdc = endTdc;
		_root = true;

		std::vector<TreeDecomposition::vertex_descriptor>::iterator end_it,
				end_it_before; // We gonna end it!

		end_it = root_vertices.begin();
		for (int i = 0; i < (root_vertices.size() - 1); i++) {
			end_it_before = end_it;
			end_it++;
//			std::cout << "ROOT: "<<*end_it_before<<" "<< *end_it << std::endl;
			do_forget(*end_it_before, *end_it);
//			std::cout << "ROOT DONE: "<<*end_it_before<<" "<< *end_it << std::endl;
		}

		/*
		 if(_tab[*end_it]->find(0)->second > _tab[*end_it]->find(1)->second)
		 solution=_tab[*end_it]->find(1)->second;
		 std::cout << "SOLUTION: " << solution << std::endl;
		 */

		BOOST_FOREACH(tabIterator i, *_tab[*end_it]) {
			std::cout << "solution array: " << i.first << " value: " << i.second
					<< std::endl;
			if (i.second < solution) {
				if (i.first == 0)
					solution = i.second;
				else
					solution = i.second + 1;
			}
		}
	}

#if PEEK_BAG
	std::cout << "Max difference was: "<< maxDiff << " with bagsize: "<< bagSize <<std::endl;
#endif

	std::cout << "SOLUTION: " << std::endl;
	std::cout << solution;

	/*
	 int tmp7 = 0;
	 tabElement test0 = forgetPosition(6,tmp7);
	 std::cout << "forgetting 110 on pos 0 " << test0;
	 //	free(_tab[t]);
	 */

}

void VertexCoverLazy::do_introduce(
		const TreeDecomposition::vertex_descriptor& child,
		const TreeDecomposition::vertex_descriptor& t) {

//	std::cout << "INTRODUCING: " << t << " INTO CHILD: " << child << std::endl;
//	std::cout << "ODD ONE OUT: " << odd_position(child, t) << std::endl;

#if OUT
	std::cout << "I";
#endif

#if DEBUG_OUT
	std::cout <<"Inserting into " << t<< " child " << child <<std::endl;
#endif

	std::vector<int>::iterator positionT;
	std::vector<int>::reverse_iterator inBag_it;
	const Bag *tee;
	tee = treedecomposition()->bag(t);

	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!

	_tabSize[t] = _tabSize[child];

	aTable *tab = new aTable;
	int position = odd_position(child, t);
	tabElement offset = _myone << position;
//	aTable* jTab = _tab[child];

	tabElement tmp;
	int debug_counter = 0;

	/*
	 for(aTable::iterator i = _tab[child]->begin();i != _tab[child]->end(); i++){
	 tmp = i->first;
	 tmp = tmp+((tmp>>position)<<position);
	 tab->insert(std::make_pair(tmp,i->second));
	 debug_counter++;
	 //		std::cout << "new TabElement after I: " << binary(i->first) << std::endl;
	 }
	 */
	tab->rehash(_tabSize[t]);

	BOOST_FOREACH(tabIterator i, *_tab[child]) {
		tmp = i.first;

//		i.first = tmp+((tmp>>position)<<position);
		tabElement x = ((_myone << (position)) - _myone);
		tabElement y = ((tmp & (~x)) << 1);
		tabElement z = y | (tmp & x);

		tab->insert(std::make_pair(z, i.second));
		debug_counter++;
#if DEBUG_OUT
		std::cout << "new TabElement after I: " << binary(i.first) << " with value "<< i.second << std::endl;
#endif
	}

#if OUT
	std::cout << "(" << debug_counter << ")";
#endif
	_tabSize[t] = debug_counter;

//	std::cout << "Memory freed in Introduce for bag: " << child << std::endl;
	delete _tab[child];
//	std::cout << "Memory freed in Introduce for list: " << child << std::endl;
	_tab[t] = tab;

	/*
	 for (int i = 0;i<_tabSize[t];i++) {
	 std::cout << binary(i) << " has value " << _tab[t][i] << std::endl;
	 }
	 */

}

void VertexCoverLazy::do_forget(
		const TreeDecomposition::vertex_descriptor& child,
		const TreeDecomposition::vertex_descriptor& t) {

#if OUT
	std::cout << "F";
#endif
#if DEBUG_OUT
	std::cout << "FORGETTING: " << t << " FROM CHILD: " << child << std::endl;
#endif
	/*
	 * this is from old introduce and gets all the vertices from and to the ODD vertex
	 * in regards to all of the elements in the bag of the CHILD!
	 *
	 * => we need to expand the table entries in the child
	 */

	std::vector<int>::iterator positionT;
	std::vector<int>::reverse_iterator inBag_it;
	const Bag *tee;
	tee = treedecomposition()->bag(child);

	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!

	unsigned long position = odd_position(child, t);
	tabElement offset = _myone << position;

	int countEdges = 0;
	std::vector<tabElement> edge_list;
//	unsigned  long* edge_list = (unsigned  long*) malloc( ((unsigned  long) tBag.size()) * sizeof(tabElement));
	bool inBag;
	int targetPos, sourcePos;

	LabeledGraph::out_edge_iterator oe_it, oe_end;

	// FIXME: Slow, maybe use set or tree?
	for (positionT = tBag.begin(); positionT < tBag.end(); positionT++) {

		//getting the out edges of all vertices in Bag
		for (tie(oe_it, oe_end) = graph()->out_edges(*positionT);
				oe_it != oe_end; oe_it++) {

			inBag = false;
			inBag_it = tBag.rbegin();
			targetPos = 0;
			// testing whether the target is in bag and which position it has
			while (inBag_it < tBag.rend()) {
				if (*inBag_it == graph()->target(*oe_it)) {
					inBag = true;
					break;
				}
				inBag_it++;
				targetPos++;
			}

			inBag_it = tBag.rbegin();

			sourcePos = 0;
			// looking up position of source in the bag
			while (inBag_it < tBag.rend()) {
				if (*inBag_it == graph()->source(*oe_it)) {
					break;
				}
				inBag_it++;
				sourcePos++;
			}
#if DEBUG_OUT_2
			std::cout << "target: " << targetPos << " source is: " << sourcePos << "position is " << position << "inbag: " << inBag << std::endl;
#endif

			if (inBag && (sourcePos == position || targetPos == position)) {
				edge_list.push_back(
						(_myone << sourcePos) | (_myone << targetPos));
#if DEBUG_OUT
				std::cout << "Inserted EDGE: " << binary(edge_list[countEdges]) << "Source is: " << graph()->source(*oe_it) << std::endl;
#endif
				countEdges++;
			}

		}

	}

	/*
	 * now we have countEdges -> is the number of found edges
	 * a matrix which has 1s which stand for the vertices connected
	 * for example:
	 * 001100 -> vertices 2 and 3 are connected through an edge
	 */

	aTable *tabChild = (_tab[child]);

	aTable *tab = new aTable;

	tabElement adjacentVertices = 0;

	for (unsigned int j = 0; j < countEdges; j++) {
#if DEBUG_OUT_2
		std::cout << "selecting edge: " <<binary(edge_list[j]) << std::endl;
#endif
		adjacentVertices = adjacentVertices | edge_list[j];
	}

//	adjacentVertices=adjacentVertices&(~offset);
#if DEBUG_OUT
	std::cout << "adjacentVertices: " << binary(adjacentVertices) << std::endl;
#endif

	unsigned long newTabsize = 0;
	tab->rehash(_tabSize[child]);

	BOOST_FOREACH(tabIterator i, *tabChild) {
#if DEBUG_OUT
		std::cout << "Processing Element: " << binary(i.first) << " with value " << i.second << std::endl;
#endif
		if ((i.first & offset) == 0) {

			// Expanding the entry to: all adjacent vertices are 1
			tabElement tmp = i.first;
			tabElement tmp3 = (tmp | adjacentVertices);
			tabElement tmp2 = forgetPosition(tmp3, position);
//			std::cout << "Position: " << position << " i.first: " << tmp << "expanded val: " << tmp3 << " adjVert: " << adjacentVertices << " newkey: " << tmp2 << std::endl;
			int second = i.second;

			// Expanding the entry into a new one with: to be forgotten vertice = 1

			// Filling in the reduced elements
			insertIntoTab(tab, std::make_pair(tmp2, second), newTabsize);

			insertIntoTab(tab,
					std::make_pair(forgetPosition(tmp, position), second + 1),
					newTabsize);

#if DEBUG_OUT
			std::cout << "Inserting into new table: " << binary(tmp2) << " and " << /*numOnes(tmp2) - numOnes(tmp) +*/second << std::endl;
			std::cout << "Inserting into new table: " << binary(forgetPosition(tmp, position)) << " and " << second+1 << std::endl;
#endif

		} else {

			insertIntoTab(tab,
					std::make_pair(forgetPosition(i.first, position),
							i.second + 1), newTabsize);

#if DEBUG_OUT
			std::cout << "Inserting into new table: " << binary(forgetPosition(i.first, position)) << " and " << i.second +1<< std::endl;
#endif
		}
	}

#if OUT
	std::cout << "(" << newTabsize << ")";
#endif

#if DOM_REMOVE
	newTabsize = newTabsize - findDominated(tab, newTabsize);
	tab->rehash(newTabsize);
#endif

	// forget the child's table
	delete _tab[child];
	_tab[t] = tab;
	_tabSize[t] = newTabsize;

#if PEEK_BAG
	int max = 0;
	int min = std::numeric_limits<int>::max();
	BOOST_FOREACH(tabIterator i, *_tab[t]) {
		if (i.second < min)
		min = i.second;
		if (i.second > max)
		max = i.second;
	}
	if (maxDiff < max - min) {
		maxDiff = max - min;
		bagSize = tBag.size();
	}

#endif

}

void VertexCoverLazy::do_join(const TreeDecomposition::vertex_descriptor& left,
		const TreeDecomposition::vertex_descriptor& right,
		const TreeDecomposition::vertex_descriptor& t) {

#if DEBUG_OUT
	std::cout << "JOINING "<<t << "with" << left <<"and"<<right <<" odd one is: " << odd_position(right, t) << std::endl;
#endif

#if TIMED
	clock_t start1, end1, start2, end2;
	start1 = clock();
#endif

	aTable *tab = new aTable;
	unsigned long countSize = 0;

#if OUT
	std::cout << "[" << _tabSize[left] << ";" << _tabSize[right] << "]";
#endif

	BOOST_FOREACH(tabIterator i, *_tab[left]) {
		BOOST_FOREACH(tabIterator j, *_tab[right]) {
			insertIntoTab(tab,
					std::make_pair((i.first | j.first), i.second + j.second),
					countSize);
		}
	}

#if 0

	/*
	 * Second loop. we discard all the tabElements from right, which have a
	 * perfect match in left.
	 * Those who are still looking for their soulmates are attending yet another
	 * matchmaking party at lonelyLovers.
	 *
	 */
#if 0
	BOOST_FOREACH(tabIterator i, *_tab[right]) {
		aTable::iterator soulMate;
		tabElement ifirst = i.first;
		soulMate = _tab[left]->find(ifirst);
		if (soulMate == _tab[right]->end()) {
			// the entry from the tab has no soulmate!
			lonelyLovers2.push_back(ifirst);
		}

	}

	for (it = lonelyLovers2.begin(); it < lonelyLovers2.end(); it++) {
#endif

		BOOST_FOREACH(tabIterator it, *_tab[right]) {

			tabElement itStar= it.first;
			int ones = numOnes(itStar);
			tabElement best= std::numeric_limits<int>::max();
			// in best we save the l1 | l2 , and l1 ^ l2 + l2.sec for faster access
			BOOST_FOREACH(tabIterator i, *_tab[left]) {
				int vgl = numOnes(i.first | itStar) + i.second;
#if DEBUG_OUT_3
				int bla = *it;
				int bla2 = i.first;

				std::cout << "We compare " << binary(bla2);
				std::cout << " and " << binary(bla);
				std::cout << " OR is " << binary(i.first | *it);
				std::cout<< " number of ones " << numOnes(i.first | *it) << " our value is " << vgl << std::endl;
#endif
				if ( vgl < best) {
					best = vgl;
				}
			}
			//best = best + numOnes(*it);

			BOOST_FOREACH(tabIterator i, *_tab[left]) {
				int vgl = ones + i.second;
				if ( vgl <= best) {
					//std::cout << ".";
					insertIntoTab(tab, std::make_pair((i.first | itStar), i.second + it.second), countSize);
				}
			}

#if DEBUG_OUT
			std::cout << "Nr. " << countSize << " was tried to be inserted with " << binary((best.first | *it)) << " and " << (best.second - numOnes(best.first ^ *it) + _tab[right]->find(*it)->second) << std::endl;
#endif

		}
#endif

#if OUT
	std::cout << "(" << countSize << ")";
#endif

#if TIMED
	end1 = clock();
	start2 = clock();
#endif

#if DOM_REMOVE
	countSize = countSize - findDominated(tab, countSize);
#endif

#if TIMED
	end2 = clock();

	std::cout << std::endl;
	std::cout << "Inserting:" << end1 << " " << start1 << std::endl;
	std::cout << "Cleaning (DOM):" << end2 << " " << start2 << std::endl;
#endif

	_tab[t] = tab;
	_tabSize[t] = countSize;
	delete _tab[right];
	delete _tab[left];

#if PEEK_BAG

	const Bag *tee;
	tee = treedecomposition()->bag(t);

	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!

	int max = 0;
	int min = std::numeric_limits<int>::max();
	BOOST_FOREACH(tabIterator i, *_tab[t]) {
		if (i.second < min)
		min = i.second;
		if (i.second > max)
		max = i.second;
	}
	if (maxDiff < max - min) {
		maxDiff = max - min;
		bagSize = tBag.size();
	}

#endif

}

void VertexCoverLazy::post_solve() {
	std::cout << std::endl;
}

std::vector<int> VertexCoverLazy::getElements(const Bag* b) {
	std::vector<int> tBag;
	tBag.reserve(b->width());
	for (int i = 0; i < b->width(); i++) {
		tBag.push_back(b->get(i));
	}
	return tBag;
}

/*
 * Takes two bags and returns one item, which is in only one of the two bags.
 *
 * -> returns -1337 if there is no difference in the two bags.
 * -> uses reverse iterator, so best if the item we are looking for is at the end of the vector.
 * -> there should be only one different item if we have a "nice" tdc.
 */
unsigned int VertexCoverLazy::odd_one_out(
		const TreeDecomposition::vertex_descriptor& child,
		const TreeDecomposition::vertex_descriptor& t) {

	std::vector<int>::reverse_iterator positionT, positionC;
	const Bag *tee, *cee;
	tee = treedecomposition()->bag(t);
	cee = treedecomposition()->bag(child);

	std::vector<int> childBag;
	childBag = getElements(cee);
	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!

	positionT = tBag.rbegin();
	positionC = childBag.rbegin();

	while (positionC < childBag.rend() || positionT < tBag.rend()) {

		if (positionC == childBag.rend())
			return *positionT;
		else if (positionT == tBag.rend())
			return *positionC;
		else if (*positionC > *positionT) {
			return *positionC;
		} else if (*positionC < *positionT) {
			return *positionT;
		} else {
			positionC++;
			positionT++;
		}
	}

	// this shouldn't happen, no difference in the bags.
	return -1337;
}

/*
 * Takes two bags and returns the positon of one item, which is in only one of the two bags.
 *
 * -> returns -1337 if there is no difference in the two bags.
 * -> there should be only one different item if we have a "nice" tdc.
 */
unsigned long VertexCoverLazy::odd_position(
		const TreeDecomposition::vertex_descriptor& child,
		const TreeDecomposition::vertex_descriptor& t) {

	unsigned long odd_position = 0;
	std::vector<int>::reverse_iterator positionT, positionC;
	const Bag *tee, *cee;
	tee = treedecomposition()->bag(t);
	cee = treedecomposition()->bag(child);

	std::vector<int> childBag = getElements(cee);
	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!

	positionT = tBag.rbegin();
	positionC = childBag.rbegin();

	while (positionC < childBag.rend() || positionT < tBag.rend()) {

		if (positionC == childBag.rend())
			return (odd_position);
		else if (positionT == tBag.rend())
			return (odd_position);
		else if (*positionC > *positionT) {
			return (odd_position);
		} else if (*positionC < *positionT) {
			return (odd_position);
		} else {
			positionC++;
			positionT++;
			odd_position++;
		}
	}

	// this shouldn't happen, no difference in the bags.
	return -1337;
}

/*
 int VertexCoverLazy::numOnes(tabElement in){
 tabElement x = in; // x has the number of '1's in in
 x = (x & 0x5555555555555555LL) + ((x >> 1) & 0x5555555555555555LL);
 x = (x & 0x3333333333333333LL) + ((x >> 2) & 0x3333333333333333LL);
 x = (x & 0x0F0F0F0F0F0F0F0FLL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FLL);
 x = (x & 0x00FF00FF00FF00FFLL) + ((x >> 8) & 0x00FF00FF00FF00FFLL);
 x = (x & 0x0000FFFF0000FFFFLL) + ((x >>16) & 0x0000FFFF0000FFFFLL);
 x = (x & 0x00000000FFFFFFFFLL) + ((x >>32) & 0x00000000FFFFFFFFLL);
 return (int)x;
 }
 */

inline unsigned long VertexCoverLazy::numOnes(const tabElement& in) {
	return mpz_popcount(in.get_mpz_t());
}

/*
 * Prints all dominated configurations to stdout
 * and deletes them.
 * slowly.
 * very slowly.
 */
int VertexCoverLazy::findDominated(aTable* tab, unsigned long projectedSize) {

	int i = 0; // counter for dominated entries.
	int diff = 0;

	boost::unordered_set<tabElement>* deletemap = new boost::unordered_set<
			tabElement>;
	deletemap->rehash(projectedSize);

	BOOST_FOREACH(tabIterator it1, *tab) {

		BOOST_FOREACH(tabIterator it2, *tab) {
			if (it1.second + numOnes(it2.first & (~it1.first)) <= it2.second
					&& (it1.first != it2.first)) {
				deletemap->insert(it2.first);
			}
		}
	}

	BOOST_FOREACH(tabElement bla, *deletemap) {

		tab->erase(bla);
		i++;

	}

#if DOM_OUT
	std::cout << "<" << i << ">";
#endif

	return i;
}

#if DOM_INSERT == 0

void VertexCoverLazy::insertIntoTab(aTable* tab, tabIterator tabEl, tabElement & tabSize) {
	tabElement first = tabEl.first;
	aTable::iterator found = tab->find(first);
	if (found == tab->end()) {
		tab->insert(tabEl);
		tabSize++;
	} else if (found->second > tabEl.second) {
		tab->erase(first);
		tab->insert(tabEl);
	}
}

#endif

#if DOM_INSERT == 1

void VertexCoverLazy::insertIntoTab(aTable* tab, tabIterator tabEl,
		unsigned long & tabSize) {

	// Check for domination, ideally should include functionality of "normal" insertIntoTab

	// We have to insert it, except if there is one entry already dominating it
	bool insert = true;
	int diff = 0;

	boost::unordered_set<tabElement>* deletemap = new boost::unordered_set<
			tabElement>;

	BOOST_FOREACH(tabIterator it2, *tab) {

#if 0
		if (tabEl.first == it2.first) {

			// Old functionality, if same element is in table we delete it and insert the new one
			if (tabEl.second < it2.second) {
				deletemap->insert(it2.first);
			}

		} else if( (tabEl.first | it2.first) == tabEl.first ) {

			// it2 is subset of tabEl
			diff = numOnes(tabEl.first ^ it2.first);
			if (tabEl.second + diff <= it2.second) {
				// it2 is being dominated

				// now kill it2
				deletemap->insert(it2.first);
			} else if(tabEl.second >= it2.second && diff != 0) {
				// tabEl has someone dominating him
				insert = false;
			}

		} else if ((tabEl.first | it2.first) == it2.first) {
			// it1 is subset of it2
			diff = numOnes(tabEl.first ^ it2.first);
			if (it2.second >= tabEl.second && diff != 0) {
				// it2 is being dominated

				// now kill it2
				deletemap->insert(it2.first);
			} else if (it2.second + diff <= tabEl.second && diff != 0) {
				// tabEl has someone dominating him
				insert = false;
			}
		}
#endif

		if (it2.second + numOnes(it2.first & (~tabEl.first)) <= tabEl.second) {
			insert = false;
		} else if (tabEl.second + numOnes(tabEl.first & (~it2.first))
				<= it2.second) {
			deletemap->insert(it2.first);
		}
	}

	if (tab->size() == 0)
		insert = true;

	BOOST_FOREACH(tabElement bla, *deletemap) {
		tab->erase(bla);
		tabSize--;
	}

	if (insert) {
		tab->insert(tabEl);
		tabSize++;
	}

}

#endif

/*
 unsigned long VertexCoverLazy::forgetPosition(unsigned long in, int& pos){
 unsigned long x = (((in>>(pos+1))<<pos)|((in<<(64-pos))>>(64-pos)));
 return x;
 }*/

tabElement VertexCoverLazy::forgetPosition(tabElement in, unsigned long& pos) {
	tabElement x = ((_myone << (pos)) - 1);
	tabElement y = (in >> 1) & (~x);
	tabElement z = y | (in & x);
	return z;
}

} //namespace

int main(int argc, char **argv) {

	std::cout << "Initializing... " << std::endl;
	sequoia::VertexCoverLazy *app = new sequoia::VertexCoverLazy();

	std::cout << "Init Done..." << std::endl;
	app->init(argc, argv);
	app->solve();

	return 0;
}
