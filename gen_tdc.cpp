#include <mindeg.h>
#include "structures/labeled_graph.h"
#include "structures/graph.h"
#include "structures/graph_printer.h"
#include "structures/graph_factory.h"
#include "structures/treedecomposition.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>


using namespace sequoia;

void usage(char *name) {
    std::cerr << "Usage: " << name << "<input tdc>" << std::endl;
    exit(1);
}


int main(int argc, char **argv) {
    sequoia::LabeledGraph g;

    unsigned int width;
    const char *gfilename = NULL;
    const char *tfilename = NULL;

    if(argc < 2) {
        std::cerr << "Error:  No input file given" << std::endl;
        usage(argv[0]);
    }
    gfilename = argv[1];
    tfilename = argv[2];

    std::ifstream ifile(argv[1]);
    if(!ifile.good()) {
        std::cerr << "Error opening file " << argv[1] << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cerr << "Loading graph file..." << std::endl;
    try {
        sequoia::GraphFactory<LabeledGraph>::load_graph(g, gfilename);
    } catch (const std::exception &e) {
        std::cerr << "Error loading graph file: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    ifile.close();
    std::cerr << "Done. Run mindegree heuristic..." << std::endl;

    sequoia::MinHeuristic<sequoia::LabeledGraph> a(g);
    a.compute();
    sequoia::TreeDecomposition *tdc = a.get();
    if (tdc == NULL) {
        std::cerr << "Error converting file " << argv[1] << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cerr << "Done, tree decomposition has width " << tdc->width() - 1
            << ".  Make nice..." << std::endl;
    tdc->make_nice();


    typedef sequoia::GraphPrinter<sequoia::TreeDecomposition> TdcPrinter;
    try {
        TdcPrinter::write_graph(tfilename, *tdc);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
