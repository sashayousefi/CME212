#include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int,int>;
  using NodeType  = typename GraphType::node_type;
  using EdgeType  = typename GraphType::edge_type;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  Point p;

  virtual void SetUp() {
    p = Point(CME212::random(), CME212::random(), CME212::random());
  }
};

class SingleNodeFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int,int>;
  using NodeType  = typename GraphType::node_type;

  //Set up Graph and Node
  GraphType graph;
  Point p;
  NodeType node;

  virtual void SetUp() {
    p = Point(CME212::random(), CME212::random(), CME212::random());
    graph.add_node(p);
    node = graph.node(0);
  }

};


class TwoGraphFixture  : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int,int>;
  using NodeType  = typename GraphType::node_type;

  //Set up Graph and Nodes
  GraphType g, g2;
};



// test Node Size
TEST_F(GraphPointFixture, NodeSize){
  EXPECT_LE( sizeof(NodeType), 16 ) << "node size > 16 bytes";
}

// Test add node
TEST_F(GraphPointFixture, AddNode){
  
  graph.add_node(p);
  EXPECT_EQ(graph.num_nodes(),1) << "Graph does not have 1 Node";
  
  NodeType node = graph.node(0);
  EXPECT_EQ(node.position() , p) << "Node position not conserved";
  EXPECT_EQ(node.index(), 0) << "Index is not 0";
  EXPECT_EQ(node.value(), 0) << "Node value is not default";
}

// Test setting node value
TEST_F(SingleNodeFixture, SetNodeValue){
  graph.node(0).value() = 2;
  EXPECT_EQ(graph.node(0).value(),2) << "New Node value from graph is incorrect";
  EXPECT_EQ(node.value(),2) << "New Node value from proxy node is incorrect";
}

// Test node removal
TEST_F(SingleNodeFixture, RemoveNode){
  EXPECT_TRUE(graph.remove_node(node)) << "Node removal returned false";
  EXPECT_EQ(graph.num_nodes(),0) << "Graph noes not have 0 Nodes after removal";
}

// Test mltiple node removal
TEST_F(GraphPointFixture, RepeatedRemoveNode){

  //add nodes
  for (int k = 0; k < 100; ++k){
    graph.add_node(Point(CME212::random(), CME212::random(), CME212::random()), k);
  }
  //check succesful add_node
  EXPECT_EQ(graph.node(0).value(),0) << "Node value set error node (0)";
  EXPECT_EQ(graph.node(53).value(),53) << "Node value set error (53)";
  EXPECT_EQ(graph.node(99).value(),99) << "Node value set error (99)";

  //check removeal
  for (unsigned k = 0; k < 50; ++k) {
    unsigned n = unsigned(CME212::random(0, graph.num_nodes()));
    /*std::cout << n << std::endl;
    graph.node(n);
    std::cout << "the issue is with node" <<std::endl;
    std::cout << n << std::endl;*/
    graph.remove_node(graph.node(n));
  }

  EXPECT_EQ( graph.num_nodes(), 50 ) << " Did not remove 50 nodes";

  //check node indices
  bool succ = true;
  for (unsigned k = 0; succ && k < 50; ++k) {
    NodeType node = graph.node(k);
    if (node.index() != k)
      succ = false;
  }
  EXPECT_TRUE(succ);

  //test clearing
  graph.clear();
  EXPECT_EQ(graph.num_nodes(), 0) << "Graph does not have 0 Nodes after removal";
}

// Test comparisons
TEST_F(TwoGraphFixture, CompareNodes){
  //Adding 50 Nodes to Graph1 and Graph2
  for (unsigned k = 0; k < 50; ++k) {
    Point p(CME212::random(), CME212::random(), CME212::random());
    g.add_node(p);
    g2.add_node(p);
  }

  EXPECT_EQ(g.num_nodes(), 50) << "g does not have 50 nodes";
  EXPECT_EQ(g2.num_nodes(), 50) << "g2 does not have 50 nodes";

  EXPECT_EQ(g2.node(23), g2.node(23)) << "G2-G2 Node comparison == error";
  EXPECT_NE(g2.node(23), g2.node(21)) << "G2-G2 Node comparison != error";
  EXPECT_NE(g2.node(23), g.node(23)) << "G2-G1 Node comparison != error";
  EXPECT_NE(g2.node(23), g.node(21)) << "G2-G1 Node comparison != error";

  EXPECT_TRUE(g.node(23) < g.node(21) || g.node(23) > g.node(21)) << "G1-G1 Node comparison < > error";
  EXPECT_TRUE(g.node(23) < g2.node(21) || g.node(23) > g2.node(21)) << "G1-G2 Node comparison < > error";

}

