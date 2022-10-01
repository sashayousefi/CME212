/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <fstream>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "SpaceSearcher.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double c = .001;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  Point initialpos;
  NodeData() : vel(0), mass(1), initialpos(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       
  double L;    
  EdgeData() : K(100), L(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0)) {
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//position struct to get update the node's position
struct compute_position{
  compute_position(double dt): _dt(dt) {}

    void operator() (Node n){
      // Compute the t+dt position
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * _dt;
    }
  private:
    double _dt;
};

//velocity struct to update the node's velocity
template <typename F>
struct compute_velocity{
  compute_velocity(double t, double dt, F& force): \
    _t(t), _dt(dt), _force(force){}

    void operator() (Node n){
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += _force(n, _t) * (_dt/n.value().mass);
    }
  private:
    double _t;
    double _dt;
    F _force;
};

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force. Parallel method.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Constraint object defining the constraint on node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports a struct
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // update position for each node
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), compute_position(dt));
  constraint(g, t);
  // Compute the t+dt velocity
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), compute_velocity<F>(t, dt, force));
  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point spring_force = Point(0, 0, 0);
    Point grav_force = Point(0, 0, -grav); //* n.value().mass;
    Point p1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Point p2 = (*it).node2().position();
      double norm_ = norm(p1 - p2);
      spring_force += -(*it).value().K * ((p1 - p2)/norm_) \
      * (norm_ - (*it).value().L);
    }
    (void) t;
    return spring_force + grav_force;
  }
};

/*Overarching force class with the default operator. The default operator 
returns a point with zero velocity at the zero x and y coordinates. */
class Force {
  public:
    virtual Point operator()(Node n, double t) const{
      (void) t;
      (void) n;
      return Point(0, 0, 0);
    }
};

/*Gravity force class. The operator returns a point with the impact of gravity 
added to the velocity.*/
class GravityForce: public Force {
  public:
    virtual Point operator() (Node n, double t) const{
      (void) t;
      return Point(0, 0, -grav) * n.value().mass;
    }
};

/*MassSpring force class. The operator returns a point with the impact 
of the spring force added to the current object.*/
class MassSpringForce: public Force{
  public:
    virtual Point operator()(Node n, double t) const {
    Point spring_force = Point(0, 0, 0);
    Point p1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      Point p2 = (*it).node2().position();
      double norm_ = norm(p1 - p2);
      spring_force += -(*it).value().K * ((p1 - p2)/norm_) \
      * (norm_ - (*it).value().L);
      (void) t;
    }
    return spring_force;
  }
};

/*Damping force class. The operator returns a point with the impact of the 
damping force added to the current object.*/
class DampingForce: public Force{
  public:
  virtual Point operator()(Node n, double t) const {
    (void) t;
    return -c * n.value().vel;
  }
};

//functor for make_combined_force
class CombinedForce {
public:
  CombinedForce(std::vector<const Force*> forces_) : forces(forces_) {}
  Point operator()(Node n, double t) {
    Point sum = Point(0, 0, 0);
    for (auto it = forces.begin(); it != forces.end(); ++it) {
      (void) t;
      const Force* f = *it;
      sum += (*f)(n, t);
    }
    return sum;
  }
private:
  std::vector<const Force*> forces;
};

/*Combines the contribution of multiple forces and applies them to the 
current object.*/
CombinedForce make_combined_force(const Force& f1, const Force& f2,
                                  const Force& f3 = Force()) {
  std::vector<const Force*> forces {&f1, &f2, &f3};
  return CombinedForce(forces);
}

/*Overarching constraint class with the default operator. The default operator 
returns no constraints. */
class Constraint {
  public:
    virtual void operator()(GraphType& g, double t) const{
      (void) t;
      (void) g;
    }
};
/*PinConstraint class with the associated operator. The constraint keeps point:
Point(0, 0, 0) and Point(1, 0, 0) fixed to hold the cloth up. */
class PinConstraint: public Constraint{
  public:
    virtual void operator()(GraphType& g, double t) const{
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        if ((*it).value().initialpos == Point(0, 0, 0)){
          (*it).position() = Point(0, 0, 0);
        }
        if ((*it).value().initialpos == Point(1, 0, 0)){
          (*it).position() = Point(1, 0, 0);
        }
      } 
    (void) t;
  }
};
/*PlaneConstraint class with the associated operator. The constraint adds a 
plane to the figure, which causes the relevant portion of the cloth to remain 
on top of the plane. */
class PlaneConstraint: public Constraint{
  public:
  float z = -0.75;
    virtual void operator()(GraphType& g, double t) const{
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Point curr_pos = (*it).position();
        Point curr_vel = (*it).value().vel;
        double dot_prod = dot(curr_pos, Point(0, 0, 1));
        //constraint invalid
        if (dot_prod < z){
          (*it).position() = Point(curr_pos.x, curr_pos.y, z);
          (*it).value().vel = Point(curr_vel.x, curr_vel.y, 0.0);
        }
      } 
    (void) t;
  }
};

/*SphereConstraint class with the associated operator. The constraint adds a 
sphere to the figure, which causes the relevant portion of the cloth to move 
and come to rest around the specified sphere. */
class SphereConstraint: public Constraint{
  public:
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
    virtual void operator()(GraphType& g, double t) const{
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Point curr_pos = (*it).position();
        Point curr_vel = (*it).value().vel;
        Point point_dist = curr_pos - c;
        double distance = norm(point_dist);
        //constraint invalid
        if (distance < r){
          (*it).position() = (r/distance)*point_dist + c;
          Point Ri = point_dist/distance;
          (*it).value().vel = curr_vel - (dot(curr_vel, Ri))*Ri;
        }
      } 
    (void) t;
  }
};

/*TearConstraint class with the associated operator. The constraint 
adds a sphere to the figureand causes the graph to tear when it 
hits the sphere. */
class TearConstraint : public Constraint {
  public:
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
    virtual void operator()(GraphType& g, double t) const{
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Point curr_pos = (*it).position();
        Point point_dist = curr_pos - c;  
        double distance = norm(point_dist);
        //constraint invalid
        if (distance < r){
          g.remove_node(*it);
        }
      } 
    (void) t;
  }
};
/*Struct to transform a node to a point*/
struct type_to_point{
  Point operator()(const Node& n){
    return n.position();
  }
};

/*SelfCollisionConstraint class with the associated operator. The constraint 
prevents the mass_spring model from passing though itself. */
struct SelfCollisionConstraint : public Constraint{
   virtual void operator()(GraphType& g, double t) const{
     (void) t;
    struct update_velocity{
      update_velocity(SpaceSearcher<Node> &_searcher, Point _range) : searcher(_searcher), 
      range(_range) {}

      void operator()(Node n){
        const Point& center = n.position();
        double radius2 = std::numeric_limits<double>::max();
        for (auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it) { //std::accumulate?
          radius2 = std::min(radius2, normSq((*e_it).node2().position() - center));
        }
        radius2 *= 0.9; 
        Box3D smallBB = Box3D(center - range*0.01, center + range*0.01);
        for (auto small_it = searcher.begin(smallBB); small_it != searcher.end(smallBB); ++ small_it){
          Node n2 = *small_it;
          Point r = center - n2.position();
          double l2 = normSq(r);
          if (n != n2 && l2<radius2) {
            n.value().vel -= (dot(r, n.value().vel)/l2) * r;
          }
        }
      }
      SpaceSearcher<Node> searcher;
      Point range;
    };
    Box3D BigBB_sizer = Box3D(thrust::make_transform_iterator(g.node_begin(), type_to_point()), \
      thrust::make_transform_iterator(g.node_end(), type_to_point()));
    Point max_val = BigBB_sizer.max();
    Point min_val = BigBB_sizer.min();
    Point range = max_val - min_val;
    Box3D BigBB = Box3D((min_val - 0.1*range), (max_val + 0.1*range));
    SpaceSearcher<Node> big_searcher = SpaceSearcher<Node>(BigBB, g.node_begin(), g.node_end(), type_to_point());
    thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), update_velocity(big_searcher, range));
  }
};

//functor for make_combined_constraint
class CombinedConstraint {
public:
  CombinedConstraint(std::vector<const Constraint*> constraints_) : \
    constraints(constraints_) {}
  virtual void operator()(GraphType& g, double t) const{
    (void) t;
    for (auto it = constraints.begin(); it != constraints.end(); ++it) {
      const Constraint* c = *it;
      (*c)(g, t);
    }
  }
private:
  std::vector<const Constraint*> constraints;
};

/*Combines the contribution of multiple constraints and applies them to the 
current object.*/
CombinedConstraint make_combined_constraint(const Constraint& c1, \
  const Constraint& c2, const Constraint& c3 = Constraint()) {
    std::vector<const Constraint*> constraints {&c1, &c2, &c3};
    return CombinedConstraint(constraints);
}