#include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"
#include "SpaceSearcher.hpp"
#include "mass_spring.hpp"
#include "MortonCoder.hpp"


// test spread_bits
TEST(MortonCoder, Spread){
  uint32_t s = detail::spread_bits(3);
  EXPECT_EQ( s, 9) << "spread of 00011 incorrect";
}

// test Space Searcher
TEST(SpaceSearcher, BasicNodeFunc){

  // Construct an empty graph
  using GraphType = Graph<int,int>;
  using Node = typename GraphType::node_type;
  GraphType graph;

  //add a 3D grid of nodes
  int size = 3;
  double step = .2;
  for( double x = 0; x <= size; x += step){
    for( double y = 0; y <= size; y += step){
      for( double z = 0; z <= size; z += step){
        graph.add_node( Point(x,y,z) );
      }
    }
  }

  // Bounding box of all the nodes
  Box3D bigbb(Point(0,0,0), Point( size+1 ,size+1,size+1));

  // Construct the Searcher
  auto n2p = [](const Node& n) { return n.position(); };
  SpaceSearcher<Node> searcher(bigbb, graph.node_begin(), graph.node_end(), n2p);

  // Spacesearcher subset of the nodes
  Box3D bb(Point(0,0,0), Point(1.1,1.1,1.1));

  //iterate through space searcher, checking position of nodes
  int count = 0;
  for (auto it = searcher.begin(bb); it != searcher.end(bb); ++it ){
    count++;
    EXPECT_LE( norm_inf( (*it).position() ) , 1.1 ) << "node outside of box included";
  }

  //Check count
  EXPECT_EQ( count, 6*6*6 ) << "Incorrect number of nodes in space searcher, expected 216"; 

}