/*
 *
 * vertexcover.h
 *
 *  Created on: Oct 18, 2011
 *      Author: fgrossmann
 */

#include <dyn_prog_solver.h>
#include <vector>
#include <string>
#include "common.h"
#include "structures/graph.h"
#include "structures/graph_structure.h"
#include "structures/graph_structure_factory.h"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <gmpxx.h>

#include <math.h>


std::size_t hash_value(const mpz_class& K) {
	size_t w = mpz_sizeinbase(K.get_mpz_t(), 2);
	mpz_class temp = K;
	size_t h, t;
	h = 0;
	for (unsigned long i = 0; i < w; i += sizeof(size_t)) {
		t = (size_t) mpz_get_ui(temp.get_mpz_t());
		temp = temp >> sizeof(size_t);
		boost::hash_combine(h, t);
	}
	return h;

}

namespace sequoia {

typedef mpz_class tabElement;

typedef std::pair<tabElement, unsigned long> tabIterator;
typedef boost::unordered_map<tabElement, unsigned long> aTable;

class VertexCoverLazy: public DynProgSolver {

public:

	void init(int argc, char **argv);

	void graphFileName(const char* graphFileName) {
		_graphFileName = graphFileName;
	}
	const char* graphFileName() {
		return _graphFileName;
	}
	const char* treedecompositionFileName() {
		return _treedecompositionFileName;
	}
	void treedecompositionFileName(const char* treedecompositionFileName) {
		_treedecompositionFileName = treedecompositionFileName;
	}
//	    void load_graph(const char* filename);

protected:

	void do_leaf(const TreeDecomposition::vertex_descriptor& t);
	void do_root(const TreeDecomposition::vertex_descriptor& t);
	void do_introduce(const TreeDecomposition::vertex_descriptor& child,
			const TreeDecomposition::vertex_descriptor& t);
	void do_forget(const TreeDecomposition::vertex_descriptor& child,
			const TreeDecomposition::vertex_descriptor& t);
	void do_join(const TreeDecomposition::vertex_descriptor& left,
			const TreeDecomposition::vertex_descriptor& right,
			const TreeDecomposition::vertex_descriptor& t);
        void pre_solve();
	void post_solve();

private:

	unsigned int odd_one_out(const TreeDecomposition::vertex_descriptor& child,
			const TreeDecomposition::vertex_descriptor& t);
	unsigned long odd_position(
			const TreeDecomposition::vertex_descriptor& child,
			const TreeDecomposition::vertex_descriptor& t);
	unsigned long numOnes(unsigned long in);
	unsigned long numOnes(const tabElement& in);
	std::vector<int> getElements(const Bag* b);

	void insertIntoTab(aTable* tab, tabIterator tabEl, unsigned long & tabSize);
	int findDominated(aTable* tab, unsigned long projectedSize);
//		unsigned long forgetPosition(tabElement in, int& pos);

	tabElement forgetPosition(tabElement in, unsigned long& pos);

	const TreeDecomposition* treedecomposition() const {
		if (_root == false) {
			return DynProgSolver::treedecomposition();
		} else {
			// now it's action time...
			return _root_tdc;
		}
	}

	// this is for experimental peeking
	bool peeked;
	int maxDiff;
	int bagSize;

	/*
	 * returns integer as string in binary format
	 */
	char* binary(unsigned long v) {
		static char binstr[17];
		int i;

		binstr[16] = '\0';
		for (i = 0; i < 16; i++) {
			binstr[15 - i] = v & 1 ? '1' : '0';
			v = v >> 1;
		}

		return binstr;
	}

	std::string binary(tabElement m) {
		return m.get_str(2);
	}

//		const GraphStructure* _graph;

	std::vector<aTable*> _tab;

	long* _tabSize;

	const TreeDecomposition* _root_tdc;

	bool _root;
	const char* _graphFileName;
	const char* _treedecompositionFileName;

	mpz_class _myone;

};

}
