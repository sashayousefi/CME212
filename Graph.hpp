#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
//template <typename V = int>
template <typename V, typename E>

class Graph {
 private:

  struct internal_node; 
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  typedef V node_value_type;
  typedef E edge_value_type;
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  //using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  //class NodeIterator;
  struct NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodelst(), edgelst(), edgesets() {}

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
      _graph_obj = nullptr;
      _node_id = 0;
    }
    
    /** Modifyable overload for node's position */
    Point& position() {
      return _graph_obj -> nodelst.at(_node_id)._position;
    }

    /** Return this node's position. */
    const Point& position() const {
      return _graph_obj -> nodelst.at(_node_id)._position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      auto it = std::find(_graph_obj-> node_ids.begin(), _graph_obj-> \
      node_ids.end(), _node_id);
      return std::distance(_graph_obj-> node_ids.begin(), it);
    }

  
    /** Return this node's value. */
    node_value_type& value(){
      return _graph_obj-> nodelst.at(_node_id)._value;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return _graph_obj-> nodelst.at(_node_id)._value;
    }

    /** Return this node's degree. */
    size_type degree() const{
      assert(*this != Node()); //asserting that the node is valid
      return _graph_obj-> nodelst.at(_node_id).incident_edge_ids.size();
    }

    /** Return the begin iterator for the incident iterator */
    incident_iterator edge_begin() const{
      return IncidentIterator(_graph_obj, _node_id, 0);
    }

   /** Return the end iterator for the incident iterator */
    incident_iterator edge_end() const{
      return IncidentIterator(_graph_obj, _node_id, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n._node_id == _node_id) && (_graph_obj == n._graph_obj);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return _node_id < n._node_id;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Pointer back to the Graph container
    Graph* _graph_obj;
    // This nodes's unique index number
    size_type _node_id;

    /** Private Constructor */
    Node(const Graph* _graph, size_type node_id)
        : _graph_obj(const_cast<Graph*>(_graph)), _node_id(node_id) {}
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodelst.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type& val = \
    node_value_type()) {
    //inserting new node_id & internal_node pair into node map
    size_type new_id = size();
    internal_node new_node = internal_node(position, val);
    nodelst.insert(std::pair<size_type,internal_node>(new_id, new_node));
    node_ids.push_back(new_id);
    return Node(this, new_id);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n._graph_obj == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i == size()){
      return Node();
    }
    assert(i < num_nodes() && i >= 0);
    return Node(this, node_ids.at(i)); 
  }
  /** Remove an node from the graph and remove all edges incident to it.
   * @param[in] a: node to be deleted. @a a's value is not modified.
   * @pre @a a is a distinct valid node of this graph
   * @return returns boolean 0 or 1 if @a a was sucessfully removed 
   * from the graph. 
   * @post Remove @a a from nodelst and nodeids.
   * @post has_node(@a a) == false.
   * @post num_nodes() = num_nodes() - 1.
   * @post num_edges() = num_edges() - incident_edges().
   *
   * Can invalidate the node @a a which is removed. Can invalidate the 
   * incident edges which are removed. 
   *
   * Complexity: No more than O(1). Nodelst.erase calls the erase function
   * on an unordered map, which is O(1) complexity. To remove the 
   * corresponding index from the indicies vector, I use the pop and swap
   * method. This method is also O(1). The element to be removed is replaced
   * with the last element of the vector, and the last element is popped off 
   * in O(1) time.
   */

  size_type remove_node(const Node& a){
    if (nodelst.find(a._node_id) == nodelst.end()){
      return 0;
    }
    while (a.degree() != 0) {
      auto it = a.edge_begin();
      remove_edge(*it);
    }
    nodelst.erase(a._node_id);
    //implementing pop and swap for vector removal
    auto it = node_ids.begin() + a.index();
    *it = node_ids.back();
    node_ids.pop_back();
    return 1;
  }

  /** Remove an node from the graph and return a valid iterator
   * @param[in] a: @a n_it: iterator pointing to the node to be deleted.
   * @pre @a n_it is a valid iterator.
   * @return returns a valid node_iterator.
   * @post has_node(@a *n_it) == false.
   * @post num_nodes() = num_nodes() - 1.
   * @post num_edges() = num_edges() - incident_edges().
   * 
   * Can invalidate the current node iterator. Must not invalidate the iterator
   * returned.
   *
   * Complexity: No more than O(1). Calls the remove_node funcition 
   * which has complexity O(1)
   */
  node_iterator remove_node(node_iterator n_it){
    assert(n_it._graph_obj); //asserting valid iterator
    if (n_it != node_end()){
      Node to_remove = *n_it;
      remove_node(to_remove);
    }
    return node_begin();
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      _graph_obj = nullptr;
      _edge_id = 0;
    }

    /** Return a node of this Edge */
    Node node1() const{
      if (_node1_idx == -1) {
        return Node(_graph_obj, _graph_obj -> edgelst.at(_edge_id).node1_id);
      }
      else{
        return Node(_graph_obj, _node1_idx);
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (_node2_idx == -1) {
        return Node(_graph_obj, _graph_obj -> edgelst.at(_edge_id).node2_id);
      }
      else{
        return Node(_graph_obj, _node2_idx);
      }
    }

    /** Find the edge length **/
    double length(){
      return norm(node1().position() - node2().position());
    }

    /** Return this edge's value. */
    edge_value_type& value(){
      return _graph_obj-> edgelst.at(_edge_id)._value;
    }

    /** Return this edge's value. */
    const edge_value_type& value() const {
      return _graph_obj-> edgelst.at(_edge_id)._value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (_graph_obj != e._graph_obj) {return false;}
      return ((node1() == e.node1() && node2() == e.node2()) ||
      (node1() == e.node2() && node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (_graph_obj == e._graph_obj){
        return (_edge_id < e._edge_id);
      }
      else{
        if (_graph_obj < e._graph_obj) {
          return true;
        }
        else{
          return false;
        }
      }
    }

    size_type index() const {
      auto it = std::find(_graph_obj-> edge_ids.begin(), _graph_obj-> \
        edge_ids.end(), _edge_id);
      return std::distance(_graph_obj-> edge_ids.begin(), it);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* _graph_obj;
    // This edge's unique index number
    size_type _edge_id;
    //private constructor for edge
    Edge(const Graph* _graph, size_type edge_id)
        : _graph_obj(const_cast<Graph*>(_graph)), _edge_id(edge_id) {}

    int _node1_idx = -1;
    int _node2_idx = -1;

    Edge(const Graph* _graph, size_type edge_id, size_type node1_idx, \
      size_type node2_idx)
        : _graph_obj(const_cast<Graph*>(_graph)), _edge_id(edge_id), 
        _node1_idx(node1_idx), _node2_idx(node2_idx) {}  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgelst.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges() && i >= 0);
    return Edge(this, edge_ids.at(i)); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //using overload operator below for comparisons
    assert(a != Node() && b != Node());
    //compares the (a, b) set with the existing map of sets
    std::set<size_type> compare_set;
    compare_set.insert(a._node_id);
    compare_set.insert(b._node_id);
    if (edgesets.find(compare_set) != edgesets.end()){
      return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(a._graph_obj == b._graph_obj);
    assert(a._node_id != b._node_id);
    size_type new_id = num_edges();
    std::set<size_type> new_set; 
    new_set.insert(a._node_id); 
    new_set.insert(b._node_id);

    if (!has_edge(a, b)) { 
      //inserting new edge_id & internal_edge pair into edge map 
      //inserting a new (a, b) pair into the edgesets map
      internal_edge new_edge{a._node_id, b._node_id};
      edgelst.insert(std::pair<size_type,internal_edge>(new_id, new_edge));
      edgesets.insert({new_set, new_id});
      //adding values to the incident edge_lst
      nodelst.at(a._node_id).incident_edge_ids.push_back(new_id);
      nodelst.at(b._node_id).incident_edge_ids.push_back(new_id);
    }
    else{
      size_type existing_edge_id = edgesets.at(new_set);
      internal_edge* existing_edge = &edgelst.at(existing_edge_id);

      if (existing_edge->node1_id == a._node_id && existing_edge->node2_id == b._node_id) {
        return Edge(this, edgesets.at(new_set));
      }
      else {
        existing_edge->node1_id = a._node_id;
        existing_edge->node2_id = b._node_id;
        return Edge(this, edgesets.at(new_set));
      }
    }
  edge_ids.push_back(new_id);
  return Edge(this, new_id);
  }

  /** Remove an edge from the graph.
   * @param[in] a: node1 of edge to be deleted. @a a's value is not modified.
   * @param[in] b: node2 of edge to be deleted. @b b's value is not modified.
   * @pre Edge between @a node1 and @b node2 is a distinct valid edge of this graph. 
   * @pre Node @a node1 and Node @b node2 are valid nodes.
   * @return returns boolean 0 or 1 if Edge(@a a, @b b) was sucessfully removed from the graph. 
   * @post Removes Edge(@a a, @b b) from edgelst, edgesets, incident_edge_ids, and edge_ids.
   * @post has_edge(Edge(@a a, @b b)) == false
   * @post num_edges() = num_edges() - 1
   *
   * Can invalidate the edge Edge(a, b) which is removed.
   *
   * Complexity: No more than O(log(num_edges) + num(incident_edges)). Edgelst.erase calls the 
   * erase function on an unordered map, which is O(1) complexity. Edgesets.erase 
   * calls the erase function on an ordered map, which is O(log num_edges). Removing
   * the incident edge in the incident edge vector is O(num (incident edges))
   * To remove the corresponding index from the indicies vector, I use the pop and swap
   * method. This method is also O(1). The element to be removed is replaced
   * with the last element of the vector, and the last vector is popped off 
   * in O(1) time. Thus the total complexity is O(log(num_edges) + num(incident_edges)).
   */

  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a,b)){
       return 0;
    }
    std::set<size_type> curr_set; 
    curr_set.insert(a._node_id); 
    curr_set.insert(b._node_id);
    
    size_type id = edgesets.at(curr_set);
    edgesets.erase(curr_set); //erase value out of edgesets
    edgelst.erase(id); //erase value out of edgelsts
    nodelst.at(a._node_id).remove_incid_edge(id);
    nodelst.at(b._node_id).remove_incid_edge(id);

    size_type idx = Edge(this, id).index();
    auto it = edge_ids.begin() + idx;
    *it = edge_ids.back();
    edge_ids.pop_back();
    return 1;
  }
  /** Remove an edge from the graph.
   * @param[in] e: Edge e to be deleted. 
   * @pre Edge @e e distinct valid edge of this graph. 
   * @return returns boolean 0 or 1 if Edge(e) was sucessfully removed from the graph. 
   * @post Removes Edge @e e from edgelst, edgesets, incident_edge_ids, and edge_ids.
   * @post has_edge(E(@e e)) == false
   * @post num_edges() = num_edges() - 1
   *
   * Can invalidate the edge Edge(e) which is removed.
   *
   * Complexity: No more than O(log(num_edges) + num(incident_edges)). 
   * Calls remove_edge(a, b) which is O(log(num_edges) + num(incident_edges)).
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an node from the graph and return a valid iterator
   * @param[in] e_it: iterator pointing to the edge to be deleted
   * @pre @a e_it is a valid iterator
   * @return returns a valid edge_iterator.
   * @post has_edge(@a *e_it) == false
   * @post Removes @e Edge *e_it from edgelst, edgesets, incident_edge_ids, and edge_ids.
   * @post num_edges() = num_edges() - 1
   * 
   * Can invalidate the current edge iterator. Must not invalidate the iterator
   * returned.
   *
   * Complexity: No more than O(log(num_edges) + num(incident_edges)). 
   * Calls remove_edge(a, b) which is O(log(num_edges) + num(incident_edges)).
   */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge to_remove = *e_it;
    remove_edge(e_it);
    return edge_begin();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodelst.clear();
    edgelst.clear();
    edgesets.clear();
    node_ids.clear();
  }

  //
  // Node Iterator
  //
  
  //struct to convert an id to the node at that id.
  struct id2node{
    id2node(const Graph *g) : _graph_obj(const_cast<Graph*>(g)) {}
    Node operator() (size_type idx){
      return _graph_obj->node(idx);
    }
    private:
      Graph* _graph_obj;
  };

/** @struct Graph::NodeIterator
  *  @brief Iterator struct for nodes in a graph. A Random Access iterator. 
  *
  *   Users can iterate through the node list and access node properties. 
  */
  struct NodeIterator:thrust::transform_iterator<id2node, thrust::counting_iterator<size_type>, Node> {
    using super_t = thrust::transform_iterator<id2node, thrust::counting_iterator<size_type>, Node>;
    NodeIterator() {} // Default ( invalid ) constructor .
    NodeIterator(const super_t& ti) : super_t{ti} {} // KLUDGE Conversion constructor.
    private:
      friend class Graph;
      Graph* _graph_obj;
      NodeIterator(const graph_type* g, thrust::counting_iterator<size_type> counter \
        ) :  super_t{counter, id2node(g)}, _graph_obj(const_cast<Graph*>(g)) {}
    };

  /*returns the beginning iterator for NodeIterator*/
  node_iterator node_begin() const {
    thrust::counting_iterator<size_type> begin(0);
    return NodeIterator(this, begin);
  } 

  /*returns the end iterator for NodeIterator*/
  node_iterator node_end() const {
    thrust::counting_iterator<size_type> end(num_nodes());
    return NodeIterator(this, end);
  }
  
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
   class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      _graph_obj = nullptr;
      _iterator_id = 0;
      _node_id = 0;
    }

    //dereference operator for the incident iterator class
    Edge operator*() const{
        Edge e = Edge(_graph_obj, _graph_obj->nodelst.at(_node_id).incident_edge_ids.at(_iterator_id));
        if (e.node1() != Node(_graph_obj, _node_id)){
          return Edge(_graph_obj, _graph_obj->nodelst.at(_node_id).incident_edge_ids.at(_iterator_id), \
            e.node2()._node_id, e.node1()._node_id);
        }
        else{
          return e;
        }
    }

    //increment to the next iterator_id
    IncidentIterator& operator++() {
      _iterator_id += 1;
      return *this;
    }

    //equality operator within the incident_iterator class 
    bool operator==(const IncidentIterator& iterator_2) const {
      return (_graph_obj == iterator_2._graph_obj && \
      _iterator_id == iterator_2._iterator_id && \
      _node_id == iterator_2._node_id);
    }

   private:
    friend class Graph;

    Graph* _graph_obj;
    size_type _node_id;
    size_type _iterator_id;

    /*private constructor for the IncidentIterator class*/
    IncidentIterator(const Graph* _graph, size_type node_id, size_type it_val) \
      : _graph_obj(const_cast<Graph*>(_graph)), _node_id(node_id), _iterator_id(it_val) {}
  };
  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      _graph_obj = nullptr;
      _iterator_id = 0;
    }

    //dereference operator for the edge iterator
    Edge operator*() const{
      return _graph_obj->edge(_iterator_id);
    }

    //increment to the next iterator_id
    EdgeIterator& operator++() {
      _iterator_id += 1;
      return *this;
    }

    //equality operator for the edge iterator
    bool operator==(const EdgeIterator& iterator_2) const {
      return (_graph_obj == iterator_2._graph_obj && \
      _iterator_id == iterator_2._iterator_id);
    }

   private:
    friend class Graph;

    Graph* _graph_obj;
    size_type _iterator_id;

    /*private constructor for the EdgeIterator class*/
    EdgeIterator(const Graph* _graph, size_type it_val)
        : _graph_obj(const_cast<Graph*>(_graph)), _iterator_id(it_val) {}
  };

  /*returns the beginning iterator for EdgeIterator*/
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  } 

  /*returns the end iterator for EdgeIterator*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }
 private:

  //internal node definition 
  struct internal_node {

  Point _position;
  node_value_type _value;
  std::vector<size_type> incident_edge_ids;
  
  //removing incident edges
  void remove_incid_edge(size_type edge_id) {
    auto it = std::find(incident_edge_ids.begin(), incident_edge_ids.end(), \
      edge_id);
    *it = incident_edge_ids.back();
    incident_edge_ids.pop_back();
    }

  internal_node(const Point& position, \
   node_value_type val) : _position(position),
     _value(val) {};

  };

  //internal edge definition
  struct internal_edge{
  size_type node1_id;
  size_type node2_id;
  edge_value_type _value;
  
  internal_edge(size_type node1, size_type node2) : node1_id(node1), \
    node2_id(node2) {};
  };

  std::unordered_map<size_type, internal_node> nodelst; //id, internal pair
  std::unordered_map<size_type, internal_edge> edgelst; //id, internal pair
  std::map<std::set<size_type>, size_type> edgesets; //set, id pair
  std::vector<size_type> node_ids; //holds the index
  std::vector<size_type> edge_ids; 
};

#endif // CME212_GRAPH_HPP