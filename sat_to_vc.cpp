/*
 * File:   sat_to_vc.cpp
 * Author: xony
 *
 * Created on August 24, 2012, 1:42 PM
 */

#include "common.h"
#include "structures/graph.h"
#include "structures/graph_structure.h"
#include "structures/graph_structure_factory.h"
#include "structures/graph_printer.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>



using namespace sequoia;

typedef std::vector<long> clauseContainer;
typedef sequoia::Graph<> myGraph;
typedef myGraph::vertex_descriptor vertex_descriptor;
typedef std::vector<vertex_descriptor> vertex_range;

/*
 * Returns a vector with vertices in the path which is ordered.
 */
vertex_range create_path(myGraph &g,long length){
    vertex_range range;
    range.push_back(g.add_vertex());
    for(int i=1; i < length; i++){
        range.push_back(g.add_vertex());
        g.add_edge(range[i-1], range[i]);
    }
    return range;
}

/*
 * Returns the vertices which interface with the main graph
 * ie the vertices of the literals
 */
vertex_range create_gadget(myGraph &g, clauseContainer clause){
    vertex_range path1, path2, path3;
    
    path1 = create_path(g, clause.size() + 2);
    path2 = create_path(g, clause.size());
    
    
    // these are the literal vertices
    for(int i = 0; i < clause.size(); i++){
        path3.push_back(g.add_vertex());
    }
    
    // connect in a triangle
    for(int i = 0; i < clause.size(); i++){
        g.add_edge(path1[i+1], path2[i]);
        g.add_edge(path2[i], path3[i]);
        g.add_edge(path1[i+1], path3[i]);
    }
    
    return path3;
}

void load_dimacs(myGraph &g, std::istream& in) {
    long numVariables;
    long numClauses;
    long clauseCount = 0;
    long currentClauseNumber = 0;
    std::vector<vertex_range> paths;
    
    while (!in.eof()) {
        char buffer[256];
        in.getline(buffer, 256);
        std::string s(buffer);
        if (s.size() == 0) continue;
        std::vector<std::string> v;
        boost::split(v, buffer, boost::is_any_of(" \t\n\r"));
        if (v[0] == "c")
            continue;
        else if (v[0] == "p") {
            //v[1] == "cnf"
            numVariables = boost::lexical_cast<long> (v[2].c_str());
            numVariables++;
            numClauses   = boost::lexical_cast<long> (v[3].c_str());
            numClauses++;
            currentClauseNumber = 0;
            
            std::cout << "found these many Variables: " << numVariables << " and these many Clauses: " << numClauses << std::endl;
            
            
            // now create paths according to number of clauses and variables
            // we need numVariables many paths which are numVariables * numClauses long
            
            for (int i = 0; i < numVariables; i++){
                paths.push_back(create_path(g, (numVariables + 1) * numClauses * 2));
            }
            
        } else { 
            // we expect numbers here with which literals are in the clause
            // remember that the line has to end with 0
         
            clauseContainer clause;
            currentClauseNumber++;
            
            std::cout << "handling clause: " << currentClauseNumber << std::endl;
             
            long pos = 0;
            long temp = boost::lexical_cast<long> (v[pos].c_str());
            while(temp != 0){
                clause.push_back(temp);
                std::cout << "inserted into clause: " << temp << std::endl;
                pos++;
                temp = boost::lexical_cast<long> (v[pos].c_str());
            }
            
            if (clause.size() % 2 != 0){
                clause.push_back(numVariables);
                std::cout << "inserted buffer into clause: " << (numVariables) << std::endl;
            }
            
            clauseCount += clause.size();
            
            // also, remember that we need n gadgets of this type
            for (int i = 0; i < numVariables + 1; i++){
                // create a gadget and connect it with the paths
                vertex_range currentClause = create_gadget(g, clause);
                
                // connect it with it's paths
                for(int j = 0; j < clause.size(); j++){
                    if(clause[j] > 0) {
                        // the variable was positive in the clause
                        g.add_edge(paths[clause[j]-1][(currentClauseNumber*2-1) + (numClauses * 2 * i)], currentClause[j]);
                        
                    } else {
                        // the variable was negative in the clause
                        g.add_edge(paths[(-clause[j]) -1][(currentClauseNumber*2-2) + (numClauses * 2 * i)], currentClause[j]);
                    }
                        
                }
            } 

        }
        
    }
    
    clauseContainer clause;
    clause.push_back(-numVariables);
    clause.push_back(-numVariables);
    
    
    clauseCount += clause.size();
    currentClauseNumber++;
            
    // also, remember that we need n gadgets of this type
    for (int i = 0; i < numVariables + 1; i++){
        // create a gadget and connect it with the paths
        vertex_range currentClause = create_gadget(g, clause);

        // connect it with it's paths
        for(int j = 0; j < clause.size(); j++){
            if(clause[j] > 0) {
                // the variable was positive in the clause
                g.add_edge(paths[clause[j]-1][(currentClauseNumber*2-1) + (numClauses * 2 * i)], currentClause[j]);

            } else {
                // the variable was negative in the clause
                g.add_edge(paths[(-clause[j]) -1][(currentClauseNumber*2-2) + (numClauses * 2 * i)], currentClause[j]);
            }

        }
    }
    
    long numVertices = g.num_vertices();
    long solution = numVertices - (( numVariables * currentClauseNumber + clauseCount + 2 * currentClauseNumber) * (numVariables + 1));
    std::cout << "ClauseCount:" << clauseCount << std::endl;
    std::cout << "ClauseNumber:" << currentClauseNumber << std::endl;
    std::cout << "Vertices:" << numVertices << std::endl;
    std::cout << "Need Vertex Cover of Size:" << solution << std::endl;
}



/*
 *
 */
int main(int argc, char** argv) {
    const char *filename = NULL;
    const char *outfilename = NULL;
    filename = argv[1];
    outfilename = argv[2];
    boost::filesystem::path p(filename);
    if (!boost::filesystem::exists(p)) {
	std::cout << "No such file or directory." << std::endl;
        return 1;
    }
    std::fstream in(filename, std::fstream::in);
    myGraph g;
    load_dimacs(g, in);
    
    
    typedef sequoia::GraphPrinter<myGraph> GraphPrinter;
    try {
        GraphPrinter::write_graph(outfilename, g);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return 0;
}