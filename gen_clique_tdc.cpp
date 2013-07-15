#include "structures/graph.h"
#include "structures/graph_printer.h"
#include "structures/treedecomposition.h"

#include <iostream>
#include <fstream>

using namespace sequoia;

int main(int argc, char **argv) {
    typedef sequoia::Graph<> Graph;

    unsigned int width;
    const char *gfilename = NULL;
    const char *tfilename = NULL;
    float edgeprob = 1.0;
    unsigned int seed = 1;

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                << " graphfile tdcfile width"
                << std::endl;
        return 1;
    }
    gfilename = argv[1];
    tfilename = argv[2];
    width = atoi(argv[3]);

    std::fstream gfile(gfilename, std::ios::out);
    if (!gfile.good()) {
        std::cerr << "Error opening output file: " << gfilename << std::endl;
        return 1;
    }
    std::fstream tfile(tfilename, std::ios::out);
    if (!tfile.good()) {
        std::cerr << "Error opening output file: " << tfilename << std::endl;
        return 1;
    }

    if (argc > 4)
        edgeprob = atof(argv[5]);
    if (argc > 5)
        seed = atoi(argv[6]);

    srand(seed);


    int numberNodes = width;

    Graph grid(numberNodes);
    //std::cout << "Probs: " << edgeprob << std::endl;

    /* due to the rand() % 1000, the probability
     * is never 1.0, but might be 0.0.
     * We therefore swap the edgeprob and the test,
     * such that edgeprob is the prob of skipping an edge.
     */
    edgeprob = 1.0 - edgeprob;

    for (int i = 0; i < numberNodes; i++) {
        for(int j = i+1; j < numberNodes; j++){
            int rval = rand() % 1000;
            if(rval/1000.0 >= edgeprob)
                grid.add_edge(i, j);
        }
    }

    sequoia::TreeDecomposition tdc;
    sequoia::TreeDecomposition::vertex_descriptor old = 0, current;

    for (int i = 0; i < numberNodes; i++) {
        current = tdc.add_vertex();
	for (int j = 0; j <= i; j++)
        	tdc.bag(current)->add(j);
        if(i > 0)
            tdc.add_edge(current, old);
        else
            tdc.root(current);
        old = current;

    }

    typedef sequoia::GraphPrinter<Graph> GraphPrinter;
    typedef sequoia::GraphPrinter<sequoia::TreeDecomposition> TdcPrinter;
    try {
        GraphPrinter::write_graph(gfilename, grid);
        TdcPrinter::write_graph(tfilename, tdc);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
