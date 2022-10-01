/**
 * @file GraphSymmetricMatrix.hpp
 * Implimentation file for treating the Graph as a MTL Matrix
 */
 // HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include <fstream>
#include "CME212/Color.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
 

// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
struct NodeBoundaryInfo{
  int is_boundary;
  float g_value;
  
  NodeBoundaryInfo() : is_boundary(0), g_value(-1) {}
};


using GraphType = Graph<NodeBoundaryInfo, double>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;

// Define a GraphSymmetricMatrix that interfaces with MTL
class GraphSymmetricMatrix{
/* * Helper function to perform multiplication. Allows for delayed
* evaluation of results.
* Assign::apply (a, b) resolves to an assignment operation such as
* a += b, a -= b, or a = b.
* @pre @a size (s) == graph.size()
* @pre @a valid(graph g) == true */

    public:
    int s;
    GraphType* g;
    //constructor for GraphSymmetricMatrix
    GraphSymmetricMatrix(int _s, const GraphType* _g) : s(_s), \
                             g(const_cast<GraphType*>(_g)) {} 
    template <typename VectorIn, typename VectorOut, typename Assign>

    void mult(const VectorIn& v, VectorOut& w, Assign) const{
      for (auto it = g->node_begin(); it != g->node_end(); ++it) {
        NodeType n = *it;
        if (n.value().is_boundary){
          Assign::apply(w[n.index()], v[n.index()]);
        }
        else{
        double running_sum = -(double)n.degree() * v[n.index()];
        for (auto adj = n.edge_begin(); adj != n.edge_end(); ++adj){
          NodeType adj_node = (*adj).node2();
          if (!adj_node.value().is_boundary){
            running_sum += v[adj_node.index()];
            }
          }
          Assign::apply(w[n.index()], running_sum);
        }
      }
    }
    /* * Matvec forwards to MTL’s lazy mat_cvec_multiplier operator */
    template <typename Vector>
    mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
    operator*(const Vector& v) const {
        return {*this, v};
    }
    private:
    // Empty!
};
/* * The number of elements in the matrix . */
inline std::size_t size(const GraphSymmetricMatrix& A){
    return A.s * A.s;
}

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
    return A.s;
}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
    return A.s;
}

/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    if (bb.contains((*it).position())) {
      g.remove_node(*it);
    }
  }
}

// Define NodeColor and NodePosition functors

/** Position functor applied to all nodes in graph @a to set u[i] to the third
 * dimensional component
 * @param[in,out] u  The current iteration of our solution vector
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        apply the position functor to achieve nodes of the form 
 *         (x_i, y_i, u[i]).
 */
struct PositionFunctor{
  mtl::dense_vector<double>* u;
  PositionFunctor(mtl::dense_vector<double>* _u) : u(_u) {};
  Point operator() (NodeType node) {
    return Point(node.position().x, node.position().y, (*u)[node.index()]);
  }
};  

/** Color functor applied to all nodes in graph @a to color the graph
 * based on the node's value and the normalization function.
 * @param[in,out] u  The current iteration of our solution vector
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        apply the color functor to color the graph.     
 */
struct ColorFunctor{
  mtl::dense_vector<double>* u;
  ColorFunctor(mtl::dense_vector<double>* _u) : u(_u) {}
  CME212::Color operator() (NodeType node) {

  //find the absolute maximum
  double max = *std::max_element((*u).begin(), (*u).end());
  double min = *std::min_element((*u).begin(), (*u).end());
  double abs_max = 0.0;
  if (abs(max) >= abs(min)){
    abs_max = abs(max);
  }
  else{
    abs_max = abs(min);
  }
  double normalized = 0.0; 
  //set normalized value to 0 to avoid division by zero error
  if (abs_max == 0.0) {
    normalized  = 0.0;
  }
  else{
  //set normalized value to the abs_value/absolute_max --> 0 <= normalized <= 1
    normalized = abs((*u)[node.index()])/abs_max;
  }
  return CME212::Color::make_heat(normalized);
  }
};  
