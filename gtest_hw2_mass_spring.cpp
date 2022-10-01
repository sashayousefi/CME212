 #include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "Graph.hpp"
#include "mass_spring.hpp"

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;  // after adding value attribute to Edge
//using GraphType = Graph<NodeData>;        // before adding value attribute to Edge

// Define Node and Edge Type
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


TEST( MassSpringFunctions, GravityExistsProblem1 ){

  // Construct an empty graph
  GraphType graph;

  //Create a point
  Point p(1);

  //add to the graph
  Node n = graph.add_node(p);

  EXPECT_EQ(p, n.position()) << "inital position incorrect";

  //take a symp euler step
  double dt = 0.001;
  double t = 0;
  for (int i = 0; i < 10; i++){
    symp_euler_step(graph, t, dt, Problem1Force());
  }
  //compare the new point to old point
  EXPECT_NE(p, graph.node(0).position()) << "position from graph did not change";
  EXPECT_NE(p, n.position()) << "position from graph did not change";
  EXPECT_EQ(n.position(), graph.node(0).position()) << "position from node and graph not equal";
}


TEST( MassSpringFunctions, GravityExistsGravityForce ){

  // Construct an empty graph
  GraphType graph;

  //Create a point
  Point p(1);

  //add to the graph
  Node n = graph.add_node(p);

  EXPECT_EQ(p, n.position()) << "inital position incorrect";

  //take a symp euler step
  double dt = 0.001;
  double t = 0;
  for (int i = 0; i < 10; i++)
    symp_euler_step(graph, t, dt, GravityForce());

  //compare the new point to old point
  EXPECT_NE(p, graph.node(0).position())<< "position from graph did not change";
  EXPECT_NE(p, n.position())<< "position from graph did not change";
  EXPECT_EQ(n.position(), graph.node(0).position())<< "position from node and graph not equal";
}

TEST( MassSpringFunctions, MakeCombinedForcePairs ){

  // Construct an empty graph
  GraphType graph;

  //Create a point
  Point p(1);

  //Add to graph
  graph.add_node(p);

  //take a symp euler step
  double dt = 0.001;
  double t = 0;

  symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),  MassSpringForce() ));
  symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),  DampingForce() ));
  symp_euler_step(graph, t, dt, make_combined_force(MassSpringForce(),  DampingForce() ));
  
  //compare the new point to old point
  EXPECT_NE(Point(1), graph.node(0).position()) << "position did not change";

}

