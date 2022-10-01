#include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "Graph.hpp"
#include "shortest_path.hpp"
#include "subgraph.hpp"


class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int, int>;
  using NodeType  = typename GraphType::node_type;
  using NodeIter = typename GraphType::node_iterator;
  using EdgeType  = typename GraphType::edge_type;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  virtual void SetUp() {
    for(int i = 0; i < 10; i++)
      points.push_back(Point(i));
  }
  
};

// Test adding node with default value
TEST_F(GraphPointFixture, DefaultNodeVal){
  graph.add_node(points[0]);
  EXPECT_EQ( graph.node(0).value(), 0 ) << "add_node does not intalize node vale with a default 0 value";
}

// Test degree function
TEST_F(GraphPointFixture, Degree){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  graph.add_edge(n0, n1);

  EXPECT_EQ(n2.degree(),0)  << "n2 degree is 0";
  EXPECT_EQ(n1.degree(),1) << "n1 degree is 1";
} 

// Test degree function death
TEST_F(GraphPointFixture, Degree2){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n3 = NodeType();
  graph.add_edge(n0, n1);

  EXPECT_DEATH(n3.degree(), ".*");
} 

// Test node iterator
TEST_F(GraphPointFixture, NodeIter){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  int iter = 0;
  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    ++iter;
  }
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}

//Test nearest node
TEST_F(GraphPointFixture, ShortestPath){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  NodeIter nearest = nearest_node(graph, Point(0));
  EXPECT_EQ(*nearest, graph.node(0)) << " error finding nearest node " ;
}
// Test incident iterator
TEST_F(GraphPointFixture, IncidentIter){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n3, n0);
  
  int iter = 0;
  for(auto ei = n0.edge_begin(); ei != n0.edge_end(); ++ei){
    EdgeType e = *ei;
    NodeType n1 = e.node1();
    assert(n1 == n0);
    ++iter;
    
  }
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}

TEST_F(GraphPointFixture, Switching){
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);

  EdgeType e3 = graph.add_edge(n3, n4);
  EXPECT_EQ(e3.node1(), n3);
  EXPECT_EQ(e3.node2(), n4);
  EdgeType e4 = graph.add_edge(n4, n3);
  EXPECT_EQ(e4.node1(), n4);
  EXPECT_EQ(e4.node2(), n3);
}

TEST_F(GraphPointFixture, BothSwitchingandIncident){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  EdgeType e3_a = graph.add_edge(n0, n3);

  EXPECT_EQ(e3_a.node1(), n0);
  EXPECT_EQ(e3_a.node2(), n3);

  EdgeType e3_b = graph.add_edge(n3, n0);

  EXPECT_EQ(e3_b.node1(), n3);
  EXPECT_EQ(e3_b.node2(), n0);

  int iter = 0;
  for(auto ei = n0.edge_begin(); ei != n0.edge_end(); ++ei){
    EdgeType e = *ei;
    NodeType n1 = e.node1();
    assert(n1 == n0);
    ++iter;
  }

  EXPECT_EQ(e3_b.node1(), n3);
  EXPECT_EQ(e3_b.node2(), n0);
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}

TEST_F(GraphPointFixture, Edge_iter){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n0, n3);
  graph.add_edge(n3, n0);

  int iter = 0;
  EdgeType curr;
  for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
    EdgeType e = *ei;
    if (e != curr) {
      ++iter;

      curr = e;
    }
  }
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}

//Test shortest_len
TEST_F(GraphPointFixture, ShortestPath_Len){
  NodeType n3 = graph.add_node(points[0]);
  NodeType n4 = graph.add_node(points[1]);
  NodeType n5 = graph.add_node(points[2]);
  NodeType n6 = graph.add_node(points[3]);
  graph.add_edge(n3, n4);
  NodeIter nearest = nearest_node(graph, Point(0));
  NodeType root = *nearest;
  EXPECT_EQ( root, graph.node(0)) << " error finding nearest node " ;
  EXPECT_EQ(1, shortest_path_lengths(graph, root));
  graph.add_edge(n4, n5);
  EXPECT_EQ(2, shortest_path_lengths(graph, root));
  graph.add_edge(n5, n6);
  EXPECT_EQ(3, shortest_path_lengths(graph, root));
  graph.add_edge(n3, n6);
  EXPECT_EQ(2, shortest_path_lengths(graph, root));
  EXPECT_EQ(1, n6.value());
}
