/**
 * @file shortest_path.cpp
 * Implimentation file for using our templated Graph to determine shortest paths.
 */

#include <unordered_set>
#include <math.h>
#include <queue>
#include <vector>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int, int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;
using IncidentIter = typename GraphType::incident_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

/* comparison functor to find the nearest node*/
struct Comparison_Op {
  const Point& point;
  Comparison_Op(const Point& p) : point(p) {};

  bool operator()(NodeType a, NodeType b) {
    Point pt_a = a.position();
    Point pt_b = b.position();
    Point::value_type dist_a = 0;
    Point::value_type dist_b = 0;
    for (int i = 0; i < 3; i++){
      dist_a += pow((point[i] - pt_a[i]), 2);
      dist_b += pow((point[i] - pt_b[i]), 2);
    }
    return (sqrt(dist_a) < sqrt(dist_b));
  }
};

/*Function to find the nearest node. Instantiates a comparison operator for 
the comparison point and then iterates over a collection of nodes while using 
a comparison functor to find the nearest node. */
NodeIter nearest_node(const GraphType& g, const Point& point) {
  if (g.num_nodes() == 0) {
    return g.node_end();
  }
  else {
    Comparison_Op point_instance = Comparison_Op(point);
    return std::min_element(g.node_begin(), g.node_end(), point_instance);
  }   
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType& g, NodeType& root)
{
  int max_distance = 0;
  std::queue<NodeType> Q;
  std::unordered_set<unsigned int> visited;
  std::vector<int> distances(g.num_nodes(), 0);
  
  distances.at(root.index()) = 0;
  visited.insert(root.index());
  Q.push(root);

  while (!Q.empty()){
    NodeType curr = Q.front();
    Q.pop();
    for (IncidentIter it = curr.edge_begin(), end = curr.edge_end(); it != end; ++it){
      NodeType neighbor = (*it).node2(); 
      if (visited.find(neighbor.index()) == visited.end()){
        visited.insert(neighbor.index());
        Q.push(neighbor);
        int curr_dist = distances.at(curr.index());
        distances.at(neighbor.index()) = curr_dist + 1;
        neighbor.value() = distances.at(neighbor.index());
        if (max_distance < distances.at(neighbor.index())){
          max_distance = distances.at(neighbor.index());
        }
      }
    }
  }
  for (unsigned int i = 0; i < g.num_nodes(); i++){
    if (distances.at(i) == 0 && root.index() != i){
      g.node(i).value() = -1;
    }
  }
  return max_distance;
}



