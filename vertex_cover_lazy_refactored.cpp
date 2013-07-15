/*
 * vertexcover.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: fgrossmann
 */

#define OUT 1
#define DEBUG_OUT 0
#define DEBUG_OUT_1 0
#define DEBUG_OUT_2 0
#define DEBUG_OUT_3 0

#define DISPLAY_BAG 0 

#define DOM_INSERT 1
#define DOM_OUT 0
#define DOM_REMOVE 0

#define PEEK_BAG 0

#define TIMED 0

#include "vertex_cover_lazy_refactored.h"
#include <stdio.h>
#include <unistd.h>
#include <limits>
#include <boost/foreach.hpp>

#include <gmpxx.h>

#include <iomanip>
#include <mindeg.h>

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

	if (treedecompositionFileName() != NULL)
		load_treedecomposition(treedecompositionFileName());


#if PEEK_BAG
	peeked = false;
	bagSize = 0;
	maxDiff = 0;
        hypothesis = true;
#endif

}

void VertexCoverLazy::pre_solve() {
    DynProgSolver::pre_solve();

    /*if (graph() == NULL)
        std::cerr << "No graph loaded." << std::endl;
    if (treedecomposition() == NULL){
	    std::cout << "Generating tree decomposition..." << std::endl;
	    MinDegreeHeuristic<sequoia::LabeledGraph> min_degree(*graph());
	    min_degree.compute();
	    DynProgSolver::treedecomposition(min_degree.get());
	}
        //generate_treedecomposition();

    std::cout << "Making tree decomposition nice..." << std::endl;
    DynProgSolver::treedecomposition()->make_nice();

    int width = treedecomposition()->width();
    std::cout << "Tree decomposition has width: " << width - 1 << " (" << width << " cops)" << std::endl;
*/


    _tab.reserve(treedecomposition()->num_vertices() );
    std::cout << "RESERVED: " << treedecomposition()->num_vertices() << " MANY SLOTS FOR _TAB" << std::endl;
    _tabSize = new long[treedecomposition()->num_vertices() ];
    _root = false;
    //long tmp = sizeof(_tab) + sizeof(aTable) * _tab.capacity(); // where T is the element type
    //std::cout << "SIZE of VECTOR:" <<tmp;
    //exit(EXIT_FAILURE);
}

void VertexCoverLazy::do_leaf(const TreeDecomposition::vertex_descriptor& t) {
#if DEBUG_OUT
	std::cout << "LEAFING: " << t << std::endl;
#endif
	aTable *tab = new aTable();
	tabElement tE = 0;
	tabIterator *tabE3 = new tabIterator(std::make_pair(tE, 0));
	tab->insert(*tabE3);
	_tab[t] = tab;
	_tabSize[t] = 1;

}

void VertexCoverLazy::do_root(const TreeDecomposition::vertex_descriptor& t) {
#if DEBUG_OUT
	std::cout << "ROOT: " << t << std::endl;
#endif
	std::vector<int>::iterator position1, position2;
	const Bag *tee;
	tee = treedecomposition()->bag(t);

	std::vector<int> tBag = getElements(tee); // YES, the one from Prison Break!
	TreeDecomposition* endTdc = new TreeDecomposition(); //tBag.size());
	TreeDecomposition::vertex_descriptor oldVertex, nuVertex;
	std::vector<TreeDecomposition::vertex_descriptor> root_vertices;
	root_vertices.reserve(tBag.size());

	long solution = std::numeric_limits<long>::max();

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

	std::cout << "TW & SOLUTION: " << std::endl;
	std::cout << (treedecomposition()->width() -1)<< std::endl;
	std::cout << solution;
        
#if PEEK_BAG
        if (!hypothesis){
            std::cout << "WRONG WRONG WRONG" << std::endl;
        }
#endif
}

void VertexCoverLazy::do_introduce(
		const TreeDecomposition::vertex_descriptor& child,
		const TreeDecomposition::vertex_descriptor& t) {


//	std::cout << "INTRODUCING: " << t << " INTO CHILD: " << child << std::endl;
//	std::cout << "ODD ONE OUT: " << odd_position(child, t) << std::endl;

#if OUT && !DEBUG_OUT
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

	tabElement tmp;
	int debug_counter = 0;

	tab->rehash(_tabSize[t]);

	BOOST_FOREACH(tabIterator i, *_tab[child]) {
		tmp = i.first;

//		i.first = tmp+((tmp>>position)<<position);
		tabElement x = ((_myone << (position)) - _myone);
		tabElement y = ((tmp & (~x)) << 1);
		tabElement z = y | (tmp & x);

		tab->insert(std::make_pair(z, i.second));
		debug_counter++;
#if DEBUG_OUT_1
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

#if OUT && !DEBUG_OUT
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
#if DEBUG_OUT_1
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
#if DEBUG_OUT_1
	std::cout << "adjacentVertices: " << binary(adjacentVertices) << std::endl;
#endif

	unsigned long newTabsize = 0;
	tab->rehash(_tabSize[child]);

	BOOST_FOREACH(tabIterator i, *tabChild) {
#if DEBUG_OUT_1
		std::cout << "Processing Element: " << binary(i.first) << " with value " << i.second << std::endl;
#endif
		if ((i.first & offset) == 0) {

			// Expanding the entry to: all adjacent vertices are 1
			tabElement tmp = i.first;
			tabElement tmp3 = (tmp | adjacentVertices);
			tabElement tmp2 = forgetPosition(tmp3, position);
#if DISPLAY_BAG
			std::cout << "Position: " << position << " i.first: " << binary(tmp) << "expanded val: " << binary(tmp3) << " adjVert: " << binary(adjacentVertices) << " newkey: " << binary(tmp2) << std::endl;
#endif
			int second = i.second;

			// Expanding the entry into a new one with: to be forgotten vertice = 1

			// Filling in the reduced elements
			insertIntoTab(tab, std::make_pair(tmp2, second), newTabsize);

			insertIntoTab(tab,
					std::make_pair(forgetPosition(tmp, position), second + 1),
					newTabsize);

#if DEBUG_OUT_1
			std::cout << "Inserting into new table: " << binary(tmp2) << " and " << /*numOnes(tmp2) - numOnes(tmp) +*/second << std::endl;
			std::cout << "Inserting into new table: " << binary(forgetPosition(tmp, position)) << " and " << second+1 << std::endl;
#endif

		} else {

			insertIntoTab(tab,
					std::make_pair(forgetPosition(i.first, position),
							i.second + 1), newTabsize);

#if DEBUG_OUT_1
			std::cout << "Inserting into new table: " << binary(forgetPosition(i.first, position)) << " and " << i.second +1<< std::endl;
#endif
		}
	}

#if OUT
	std::cout << "(" << newTabsize << ")" << std::flush;
#endif

#if DOM_REMOVE
	newTabsize = newTabsize - findDominated(tab, newTabsize);
	tab->rehash(newTabsize);
#endif


#if DISPLAY_BAG
	int c=1;
	std::cout << std::endl;
	std::cout << "Child Table:" << std::endl;
	BOOST_FOREACH(tabIterator i, *_tab[child]) {
	    std::cout << std::setw(2) << std::right << c << ": " << std::setw(16) << binary(i.first) << " | " << i.second  << std::endl;
	    c++;
	}
#endif



	// forget the child's table
	delete _tab[child];
	_tab[t] = tab;
	_tabSize[t] = newTabsize;

#if DISPLAY_BAG
	c=1;
	std::cout << std::endl;
	std::cout << "Parten Table:" << std::endl;
	BOOST_FOREACH(tabIterator i, *_tab[t]) {
	    std::cout << std::setw(2) << std::right << c << ": " << std::setw(16) << binary(i.first) << " | " << i.second  << std::endl;
	    c++;
	}
#endif


#if PEEK_BAG
        int tw = 0;
	int max = 0;
	int min = std::numeric_limits<int>::max();
	BOOST_FOREACH(tabIterator i, *_tab[t]) {
		if (i.second < min)
		min = i.second;
		if (i.second > max)
		max = i.second;
	}
        
	if (max - min > (tBag.size() / 2) + 1) {
		std::cout << "HYPOTHESIS VOID" << std::endl;
		hypothesis = false;
	}
	if (newTabsize > pow(2.0,(tBag.size() / 2) )){
		hypothesis = false;
		std::cout << "Width: " << tBag.size() << " with tabSize " << newTabsize;
		exit(EXIT_FAILURE);
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
#if OUT && !DEBUG_OUT
	std::cout << "J";
#endif
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

#if 1
        /*
         * This is the old matchmaking orgy
         */
	BOOST_FOREACH(tabIterator i, *_tab[left]) {
		BOOST_FOREACH(tabIterator j, *_tab[right]) {
			insertIntoTab(tab,
					std::make_pair((i.first | j.first), i.second + j.second),
					countSize);
		}
	}
#endif


#if 0
        /*
         * This is waaaay more selective
         */

        // copy right table to a new table
        aTable *lonelyLovers = new aTable(*_tab[right]);
        std::vector<tabIterator> perfectMatches;
        boost::unordered_set<tabElement>* deletemap = new boost::unordered_set<
			tabElement>;

        unsigned long bestDistance = std::numeric_limits<unsigned long>::max();
        unsigned long newDistance;

        /* First loop
         * select perfect matches
         */
        BOOST_FOREACH(tabIterator i, *_tab[left]) {
            // find the perfect matches
            BOOST_FOREACH(tabIterator j, *_tab[right]) {
                newDistance = j.second * 2 + numOnes( /* j.first & (~i.first)*/  (i.first & j.first) ^ (i.first | j.first) );
                if (newDistance <= bestDistance ){
                    // push these onto a stack
                    perfectMatches.push_back(j);
                    bestDistance=newDistance;
                }
            }
            bestDistance = std::numeric_limits<unsigned long>::max();

            // take those from the stack and let them mate with the current entry
            tabIterator it;
            while (!perfectMatches.empty()){
                it = perfectMatches.back();
                insertIntoTab(tab, std::make_pair((i.first | it.first), i.second + it.second),
					countSize);
                deletemap->insert(it.first);
                perfectMatches.pop_back();
            }


        }

        // now delete all who have a match at lonelyLovers
        /*BOOST_FOREACH(tabElement bla, *deletemap) {
		lonelyLovers->erase(bla);
	}*/


        bestDistance = std::numeric_limits<unsigned long>::max();
        // now find matches for the remaining ones in lonelylovers
        BOOST_FOREACH(tabIterator i, *_tab[right]) {
            // find the perfect match
            BOOST_FOREACH(tabIterator j, *_tab[left]) {
                newDistance = j.second * 2  +  numOnes( /*j.first & (~i.first)*/ (i.first & j.first) ^ (i.first | j.first)  );
                if (newDistance <= bestDistance ){
                    // push these onto a stack
                    perfectMatches.push_back(j);
                    bestDistance=newDistance;
                }
            }
            bestDistance = std::numeric_limits<unsigned long>::max();

            // take those from the stack and let them mate with the current entry
            tabIterator it;
            while (!perfectMatches.empty()){
                it = perfectMatches.back();
                insertIntoTab(tab, std::make_pair((i.first | it.first), i.second + it.second),
					countSize);
                perfectMatches.pop_back();
            }

        }

	delete lonelyLovers;

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
       /* if (max - min > (tBag.size() / 2) + 1) {
            std::cout << "HYPOTHESIS VOID IN JOINS" << std::endl;
            hypothesis = false;
        }*/

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
