/*
 * vertexcover.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: fgrossmann
 */


#include "vertex_cover.h"
#include <stdio.h>
#include <unistd.h>
#include <limits>

namespace sequoia {




void VertexCover::init(int argc, char **argv) {


	std::cout << "Size is: " << sizeof(int) << std::endl;

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
				std::cerr << "Error: Unknown argument: " << argv[optind] << std::endl;
		}
	}


	if(argc - optind > 0) {
		std::cerr << "ERROR:  Not enough options specified" << std::endl;
	}



	std::cout << "Loading Graph from: " << graphFileName() << std::endl;

	load_graph(graphFileName());

//	std::cout << "Binary form of 2 is :" << binary(2) << std::endl;

	if (treedecompositionFileName() != NULL)
		load_treedecomposition(treedecompositionFileName());
	else
		generate_treedecomposition();


//    DynProgSolver::init();




/*
 *  FOR OLD FUNCTION
    std::cout << "Number of Vertices: " << _graph->graph()->num_vertices() << std::endl;
    std::cout << "Number of Edges: " << _graph->graph()->num_edges() << std::endl;
    std::cout << "Number of Labels: " << _graph->graph()->num_labels() << std::endl;


	WITH THE INHERITED ONE
    std::cout << "Number of Vertices: " << graph()->num_vertices() << std::endl;
    std::cout << "Number of Edges: " << graph()->num_edges() << std::endl;
    std::cout << "Number of Labels: " << graph()->num_labels() << std::endl;
    std::cout << "Out Edges of 1: " << graph()->out_degree(1) << std::endl;
    std::cout << "In Edges of 1: " << graph()->in_degree(1) << std::endl;

*/
}

void VertexCover::pre_solve() {
    DynProgSolver::pre_solve();

    _tab.reserve(treedecomposition()->num_vertices() ); // The first bag is empty
    _tabSize = new long[treedecomposition()->num_vertices() ] ;
}

void VertexCover::post_solve() {
	std::cout << std::endl;
}


void VertexCover::do_leaf(const TreeDecomposition::vertex_descriptor& t){

//	std::cout << "LEAFING: " << t << std::endl;

	int* pData = (int*) malloc (sizeof(int));
        pData[0] = 0;
        _tab[t] = pData;
        _tabSize[t] = 1;
        std::cout <<"("<< _tabSize[t] <<")";

}


void VertexCover::do_root(const TreeDecomposition::vertex_descriptor& t){

//	std::cout << "ROOT: " << t << std::endl;


	int temp = _tab[t][0];
//	std::cout << binary(0) << " has value " << _tab[t][0] << std::endl;

	for (unsigned int i = 1;i<_tabSize[t];i++) {
//		std::cout << binary(i) << " has value " << _tab[t][i] << std::endl;
		if (_tab[t][i]<temp){
			temp=_tab[t][i];

		}
	}



	std::cout << std::endl;
	std::cout << "SOLUTION: " << std::endl;
	std::cout << temp;

	free(_tab[t]);

}




void VertexCover::do_introduce(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t){


//	std::cout << "INTRODUCING: " << t << " INTO CHILD: " << child << std::endl;
//	std::cout << "ODD ONE OUT: " << odd_one_out(child, t) << std::endl;



	std::cout <<"I";

	std::vector<int>::iterator positionT;
	std::vector<int>::reverse_iterator inBag_it;
	const Bag *tee;
	tee = treedecomposition()->bag(t);




	std::vector<int> tBag ( getElements(tee) ); // YES, the one from Prison Break!
/*
	if (tBag.size() == 1){
		 * This will be a "true" leaf, with only one element in the bag
		 *
		 * we only need to initialize a tab with 2 entries: 0 and 1.
		 * return after that because then our work is done!
		 *
//		std::cout << "this happens once " << std::endl;
		int* pData = (int*) malloc (2*sizeof(int));

//		std::cout << "Memory allocated for Nr: " << t << std::endl;

		pData[0] = 0;
		pData[1] = 1;
		_tab[t] = pData;
		_tabSize[t] = 2;
		std::cout <<"("<< _tabSize[t] <<")";
		return;
	}
*/

	int position = odd_position(child, t);
	unsigned long long offset = 1L<<position;

	int countEdges = 0;
	unsigned long long* edge_list = (unsigned long long*) malloc( ((unsigned long long) tBag.size()) * tBag.size() * sizeof(long long));
//	std::cout << "Size of Bag: " << tBag.size() << std::endl;
	bool inBag;
	int targetPos, sourcePos;

	LabeledGraph::out_edge_iterator oe_it, oe_end;

	// FIXME: Slow, maybe use set or tree?
	for (positionT = tBag.begin(); positionT<tBag.end(); positionT++) {

		//getting the out edges of all vertices in Bag
		for (tie(oe_it, oe_end) = graph()->out_edges(*positionT); oe_it != oe_end; oe_it++) {

			inBag = false;
			inBag_it = tBag.rbegin();
			targetPos=0;
			// testing whether the target is in bag and which position it has
			while(inBag_it < tBag.rend()) {
				if (*inBag_it ==  graph()->target(*oe_it)){
					inBag = true;
					break;
				}
				inBag_it++;
				targetPos++;
			}

			inBag_it = tBag.rbegin();

			sourcePos=0;
			// looking up position of source in the bag
			while(inBag_it < tBag.rend()) {
				if (*inBag_it ==  graph()->source(*oe_it)){
					break;
				}
				inBag_it++;
				sourcePos++;
			}

			if(inBag && (((1L<<position)&(1L << sourcePos) | (1L << targetPos))!=0)) {
				edge_list[countEdges] = (1L << sourcePos) | (1L << targetPos);
//				std::cout << "Inserted: " << binary(edge_list[countEdges]) << "Source is: " << graph()->source(*oe_it) << std::endl;
				countEdges++;
			}

		}

	}

	bool validVertexCover0;
	bool validVertexCover1;
	unsigned long long j, k1, k2; //k1 and k2 are the target positions of i in the new table.
	int modCounter = 0;


	_tabSize[t]= 2*_tabSize[child];
	std::cout <<"("<< _tabSize[t] <<")";
	int* pData = (int*) malloc (((long long)_tabSize[t])*sizeof(int));

//	std::cout << "Memory allocated for Nr: " << t << std::endl;
	for (unsigned long long i=0;i<_tabSize[child];i++) {
		
		int x=(1<<position) -1;
		int y = ((i & (~x)) << 1);
		int z = y | (i & x);

		//k1=i+((i>>position)<<position);
		k1=z;
		k2=k1|offset;
		validVertexCover0 = true;
		j=0;
		while(j < countEdges) {


			if ((edge_list[j] & (k1)) == 0){
				validVertexCover0 = false;
				break;
			}
			j++;
		}

		validVertexCover1 = true;
		j=0;
		while(j < countEdges) {
//			std::cout << binary(k+1);
//			std::cout << " Debug-> " << binary(edge_list[j]) << " was compared to: " << (k+1) << " result was "<< (edge_list[j] & (k+1)) << std::endl;
			if ((edge_list[j] & (k2)) == 0){
				validVertexCover1 = false;
				break;
			}
			j++;
		}

		if( validVertexCover0 && _tab[child][i] !=  std::numeric_limits<int>::max()){
			pData[k1] = _tab[child][i];
//			std::cout << binary(k1) << " was set to: " << pData[k1] << std::endl;
		} else {
			pData[k1] = std::numeric_limits<int>::max();
//			std::cout << binary(k1) << " was set to: oo " << std::endl;
		}

		if( validVertexCover1 && _tab[child][i] !=  std::numeric_limits<int>::max()){
			pData[k2] = (_tab[child][i]+1);
//			std::cout <<  binary(k2)<< " was set to: " << pData[k2] << std::endl;
		} else {
			pData[k2] = std::numeric_limits<int>::max();
//			std::cout << binary(k2) << " was set to: oo " << std::endl;
		}


	}
//	std::cout << "Memory freed in Introduce for bag: " << child << std::endl;
	free(_tab[child]);
//	std::cout << "Memory freed in Introduce for list: " << child << std::endl;
	free(edge_list);
	_tab[t] = pData;


/*	for (int i = 0;i<_tabSize[t];i++) {
		std::cout << binary(i) << " has value " << _tab[t][i] << std::endl;
	}*/


}


void VertexCover::do_forget(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t){


	std::cout <<"F";
//	std::cout << "FORGETTING: " << t << " FROM CHILD: " << child << std::endl;

	unsigned long long i = 0; //counter for old table
	unsigned long long j = 0; //counter for new table

	unsigned long long position = 1L<<odd_position(child, t); // position of child in table from the right
//	std::cout << "Position is: " << odd_position(child, t) << std::endl;

	// we do need a new table now, but only half the size from the child's one
	_tabSize[t]= _tabSize[child]>>1;
	std::cout <<"("<< _tabSize[t] <<")";
//	std::cout << "New Tabsize: " << _tabSize[t] << std::endl;
	int* pData = (int*) malloc (((long long)_tabSize[t])*sizeof(int));
//	std::cout << "Memory allocated for Nr: " << t << std::endl;

	while(i < _tabSize[child]){
		do {
			i++;
			j++;
			pData[j-1]=_tab[child][i-1];
//			std::cout <<"SETTING: "<< j-1 << " FROM " << i-1 << " with value: " << pData[j-1] << std::endl;
		}  while ((i % position) != 0);

//		pData[j]=_tab[child][i];
//		std::cout <<"SETTING: "<< j << " FROM " << i << " with value: " << pData[j] << std::endl;

		j=j-position; //reset counter for new table
//		j++;

		do {
			i++;
			j++;
//			std::cout <<"CHECKING: "<< j-1 << " WITH " << i-1  << std::endl;
			if (pData[j-1]>_tab[child][i-1]){
				pData[j-1]=_tab[child][i-1];
//				std::cout <<"OVERWRITING: "<< j-1 << " FROM " << i-1 << " with value: " << pData[j-1] << std::endl;
			}
		} while ((i % position) != 0);



	}

	// forget the child's table
//	std::cout << "Memory freed for in forget child: " << child << std::endl;
	free(_tab[child]);
	_tab[t] = pData;


/*	for (int i = 0;i<_tabSize[t];i++) {
		std::cout << binary(i) << " has value " << _tab[t][i] << std::endl;
	}*/


}

void VertexCover::do_join(const TreeDecomposition::vertex_descriptor& left,
        const TreeDecomposition::vertex_descriptor& right,
        const TreeDecomposition::vertex_descriptor& t) {
//	std::cout << "JOINING "<<t << "with" << left  <<"and"<<right << std::endl;

	std::cout <<"J";
	_tabSize[t]= _tabSize[left];
	std::cout <<"("<< _tabSize[t] <<")";
	//	std::cout << "New Tabsize: " << _tabSize[t] << std::endl;
	int* pData = (int*) malloc (( (long long)_tabSize[t])*sizeof(int));

	for(unsigned long long i=0;i<_tabSize[t];i++){

		unsigned long long x = i; // x has the number of '1's in i
		x = (x & 0x5555555555555555LL) + ((x >> 1) & 0x5555555555555555LL);
		x = (x & 0x3333333333333333LL) + ((x >> 2) & 0x3333333333333333LL);
		x = (x & 0x0F0F0F0F0F0F0F0FLL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FLL);
		x = (x & 0x00FF00FF00FF00FFLL) + ((x >> 8) & 0x00FF00FF00FF00FFLL);
		x = (x & 0x0000FFFF0000FFFFLL) + ((x >>16) & 0x0000FFFF0000FFFFLL);
		x = (x & 0x00000000FFFFFFFFLL) + ((x >>32) & 0x00000000FFFFFFFFLL);

//		std::cout << "Binary: " << binary(i) << " has ones: " << x << std::endl;

		if (_tab[right][i] == std::numeric_limits<int>::max()|| _tab[left][i]== std::numeric_limits<int>::max())
			pData[i] =std::numeric_limits<int>::max();
		else
			pData[i] = _tab[right][i] + _tab[left][i] - x;

	}

	_tab[t] = pData;
	free(_tab[left]);
	free(_tab[right]);


}

/*
void VertexCover::load_graph(const char* filename) {
    _graph = GraphStructureFactory::load_leda(filename);
    DynProgSolver::graph(_graph->graph());
    std::cout << "Loaded Graph of order " << _graph->graph()->num_vertices() << " with Vocabulary: "
            << _graph->vocabulary()->toString() << std::endl;
}
*/

/*
 * Takes two bags and returns one item, which is in only one of the two bags.
 *
 * -> returns -1337 if there is no difference in the two bags.
 * -> uses reverse iterator, so best if the item we are looking for is at the end of the vector.
 * -> there should be only one different item if we have a "nice" tdc.
 */
int VertexCover::odd_one_out(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t){

	std::vector<int>::reverse_iterator positionT, positionC;
	const Bag *tee, *cee;
	tee = treedecomposition()->bag(t);
	cee = treedecomposition()->bag(child);

	std::vector<int> childBag ( getElements(cee) );
	std::vector<int> tBag ( getElements(tee) ); // YES, the one from Prison Break!

	positionT = tBag.rbegin();
	positionC = childBag.rbegin();

	while(positionC < childBag.rend() || positionT < tBag.rend()) {

		if (positionC == childBag.rend())
			return *positionT;
		else if (positionT == tBag.rend())
			return *positionC;
		else if (*positionC > *positionT){
			return *positionC;
		}
		else if (*positionC < *positionT ) {
			return *positionT;
		}
		else {
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
int VertexCover::odd_position(const TreeDecomposition::vertex_descriptor& child, const TreeDecomposition::vertex_descriptor& t){

	int odd_position = 0;
	std::vector<int>::reverse_iterator positionT, positionC;
	const Bag *tee, *cee;
	tee = treedecomposition()->bag(t);
	cee = treedecomposition()->bag(child);

	std::vector<int> childBag ( getElements(cee) );
	std::vector<int> tBag ( getElements(tee) ); // YES, the one from Prison Break!

	positionT = tBag.rbegin();
	positionC = childBag.rbegin();

	while(positionC < childBag.rend() || positionT < tBag.rend()) {

		if (positionC == childBag.rend())
			return (odd_position);
		else if (positionT == tBag.rend())
			return (odd_position);
		else if (*positionC > *positionT){
			return (odd_position);
		}
		else if (*positionC < *positionT ) {
			return (odd_position);
		}
		else {
			positionC++;
			positionT++;
			odd_position++;
		}
	}

	// this shouldn't happen, no difference in the bags.
	return -1337;
}

std::vector<int> VertexCover::getElements(const Bag* b) {
	std::vector<int> tBag;
	tBag.reserve(b->width());
	for (int i = 0; i < b->width(); i++) {
		tBag.push_back(b->get(i));
	}
	return tBag;
}



} //namespace

int main(int argc, char **argv) {
	std::cout << "Initializing... " << std::endl;
    sequoia::VertexCover *app = new sequoia::VertexCover();

    std::cout << "Init Done..." << std::endl;
    app->init(argc, argv);
    app->solve();
    return 0;
}


