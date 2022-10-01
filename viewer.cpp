/**
 * @file viewer.cpp
 * Test script for the SFML_Viewer and Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point list
 *
 * Prints
 * A B
 * where A = number of nodes
 *       B = number of edges
 * and launches an SFML_Viewer to visualize the system.
 */
 
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

bool check_file_exists(std::ifstream& file, std::string arg){
  if (!file.good()){
    std::cerr << arg << " provided is not a valid filename\n";
    return false;
  }
  return true;
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph<int, int>;
  //using GraphType = Graph;
  
  using NodeType  = typename GraphType::node_type;
  //using EdgeType = typename GraphType::edge_type;

  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);

  // Check that the files specified by the input arguments exist, and exit if not
  if ( !check_file_exists(nodes_file, "NODES_FILE") || !check_file_exists(tets_file, "TETS_FILE")){
    exit(1);
  }

  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);


  // Print number of nodes and edges
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;


  // Launch a viewer
  CME212::SFML_Viewer viewer;

  //viewer.draw_graph_nodes(graph);  // Draw only the nodes
  //viewer.draw_graph(graph);      // Draw the nodes and edges

  // Center the view and enter the event loop for interactivity
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  viewer.event_loop();

  return 0;
}