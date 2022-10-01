/**
 * @file poisson.cpp
 * Test script for treating for using the GraphSymmetricMatrix class
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */
#include <iostream>
#include "CME212/SFML_Viewer.hpp"
#include "GraphSymmetricMatrix.hpp" 
#include <math.h>

// Define visual_iteration that inherits from cyclic_iteration
/** Traits that MTL uses to determine properties of our GraphSymmetricMatrix.*/
namespace mtl {
namespace ashape {

/* * Define GraphSymmetricMatrix to be a non - scalar type.*/
template <>
struct ashape_aux<GraphSymmetricMatrix> {
    typedef nonscal type;
    };
} // end namespace ashape

/** GraphSymmetricMatrix implements the Collection concept
* with value_type and size_type*/
template <>
struct Collection<GraphSymmetricMatrix> {
    typedef double value_type;
    typedef unsigned size_type;
    };
} // end namespace mtl

/** @class VisualIteration
 * @brief VisualIteration class that interfaces with viewer.
 *
 * Visual Iteration contains the update_viewer method, which updates the viewer
 * at each step of the iteration. It uses the color and position functors to
 * re-color and re-draw the graph at each converging step.
 */
namespace itl{
  template <class Real, class OStream = std::ostream>
  class visual_iteration : public itl::cyclic_iteration<Real> {
    public:
    typedef cyclic_iteration<Real> super;
    typedef visual_iteration self;
    
    /*visual iteration constructor*/
    visual_iteration(const GraphType* g, CME212::SFML_Viewer* viewer, \
                    ColorFunctor func_color, PositionFunctor func_pos,\
                    const mtl::dense_vector<double>& b, int max_iter_, \
                    Real tol_, Real atol_ = Real(0), int cycle_ = 50,\
                    OStream& out = std::cout) \
                    : super(b, max_iter_, tol_, atol_, cycle_, out),\
                    graph(const_cast<GraphType*>(g)), \
                    viewer(viewer), colfunc(func_color), posfunc(func_pos) {}

    /* checks whether convergnce or max iter has been reached*/
    bool finished(){
      update_viewer();
      return super::finished();
    }

    template <typename T>
    /* checks whether convergnce or max iter has been reached*/
    bool finished(const T& r){
      update_viewer();
      bool ret = super::finished(r);
      return ret;
    }

    /*clears the graph and applies the new colors and postiions for 
    the next iteration.*/
    void update_viewer(){
      auto node_map = viewer->empty_node_map(*graph);
      viewer->clear();
      viewer->add_nodes(graph->node_begin(), graph->node_end(), colfunc, \
                        posfunc, node_map);
      viewer->add_edges(graph->edge_begin(), graph->edge_end(), node_map);
    }
    protected:
      GraphType* graph;
      CME212::SFML_Viewer* viewer;
      ColorFunctor colfunc;
      PositionFunctor posfunc;
  };
} // namespace itl

/** Traits that MTL uses to determine properties of our IdentityMatrix.*/
int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  // update the node's g_value and boolean is_boundary value.  
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    NodeType n = *it;
    Point p1 = Point(-0.6, -0.2, -1);
    Point p2 = Point(0.6, 0.2, 1);
    Box3D BB = Box3D(p1, p2);
    if (norm_inf(n.position()) == 1){
      n.value().g_value = 0;
      n.value().is_boundary = 1;
    }
    else if ((norm_inf(n.position() - Point(0.6, 0.6, 0)) < 0.2) \
    or (norm_inf(n.position() - Point(-0.6, 0.6, 0)) < 0.2) \
    or (norm_inf(n.position() - Point(0.6, -0.6, 0)) < 0.2) \
    or (norm_inf(n.position() - Point(-0.6, -0.6, 0)) < 0.2)) {

      n.value().g_value = -0.2;
      n.value().is_boundary = 1;
    }
    else if (BB.contains(n.position())){
      n.value().g_value = 1.0;
      n.value().is_boundary = 1;
    }
    else{
      n.value().g_value = -1.0;
      n.value().is_boundary = 0;
    }
  } 

  // using MTL's conjugate gradient solver
  int N = graph.size(); 

  typedef GraphSymmetricMatrix matrix_type;
  // Set up a matrix
  matrix_type A(N, &graph);

  // Create an preconditioner
  itl::pc::identity<matrix_type> L(A);

  // Set b according to boundary position. Start with x == 1.
  mtl::dense_vector<double> x(N, 1.0), b(N, 0.0);
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    NodeType n = *it;
    if (n.value().is_boundary){
      b[n.index()] += n.value().g_value;
    }
    else{
      double h2 = pow(graph.edge(0).length(), 2);
      b[n.index()] += h2*5*std::cos(norm_1(n.position())); 
      for (auto adj = n.edge_begin(); adj != n.edge_end(); ++adj){ 
        NodeType adj_node = (*adj).node2();
        if (adj_node.value().is_boundary){
          b[n.index()] -= adj_node.value().g_value;
        }
      }
    }
  }

  assert(b != x);
  // Termination criterion: r < 1e-10 * b or N iterations
  //itl::cyclic_iteration<double> iter(b, 300, pow(10, -10), pow(10, -10), 50);

  // Solve Ax == b with left preconditioner P
  //itl::cg(A, x, b, L, iter);

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  // Center the view and enter the event loop for interactivity
  //bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]() {
  itl::visual_iteration<double> vis_iter(&graph, &viewer, ColorFunctor(&x), \
                     PositionFunctor(&x), 
                     b, 300, pow(10, -10), 
                     pow(10, -10), 50);


  // Solve Ax == b with left preconditioner P
  itl::cg(A, x, b, L, vis_iter);  
  });  // simulation thread

  viewer.event_loop();

  return 0;
}
