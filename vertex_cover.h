/*
 *
 * vertexcover.h
 *
 *  Created on: Oct 18, 2011
 *      Author: fgrossmann
 */


#include "dyn_prog_solver.h"
#include <vector>
#include "common.h"
#include "structures/graph.h"
#include "structures/graph_structure.h"
#include "structures/graph_structure_factory.h"

namespace sequoia {

class VertexCover : public DynProgSolver {


	public:


		void init(int argc, char **argv);

	    void graphFileName(const char* graphFileName) { _graphFileName = graphFileName; }
	    const char* graphFileName() { return _graphFileName; }
	    const char* treedecompositionFileName() { return _treedecompositionFileName; }
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

		int odd_one_out(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t);
		int odd_position(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t);
                std::vector<int> getElements(const Bag* b);
                
		/*
		 * returns integer as string in binary format
		 */
		char* binary (int v) {
			static char binstr[17] ;
			int i ;

			binstr[16] = '\0' ;
			for (i=0; i<16; i++) {
				binstr[15-i] = v & 1 ? '1' : '0' ;
				v = v / 2 ;
			}

			return binstr ;
		}

//		const GraphStructure* _graph;

		std::vector<int*> _tab;
		long* _tabSize;
		int _hilf;


	    const char* _graphFileName;
	    const char* _treedecompositionFileName;




};

}
