/**
 * @file subgraph.hpp
 * Implimentation file for viewing a subgraph from our Graph
 */


#include <fstream>
#include <iterator>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>>
{
 public:
  // Get all of the iterator traits and make them our own
  using value_type        = typename std::iterator_traits<It>::value_type;
  using pointer           = typename std::iterator_traits<It>::pointer;
  using reference         = typename std::iterator_traits<It>::reference;
  using difference_type   = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  typedef filter_iterator<Pred, It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {}

    /* dereference operator for the filter iterator */ 
    value_type operator*() {
      return *it_;
    }

    /* ++operator for the filter iterator - increment operator. 
    skips nodes that do not satisfy a predicate*/ 
    filter_iterator& operator++(){
      do {
        ++it_;
      } while (it_ != end_ && !valid_node()); 
      return *this;
    }

    /* equality operator for the filter iterator */ 
    bool operator==(const self_type& iterator_2) const{
      return iterator_2.end_ == end_ && iterator_2.it_ == it_ \
        && typeid(p_) == typeid(iterator_2.p_);
    }

    /* checking if node is valid */ 
    bool valid_node() {
      return p_(*it_);
    }

 private:
  Pred p_;
  It it_;
  It end_;

};


/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

/** Interesting predicate for HW1 #4 */
struct InterestingPredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return norm(n.position()) < 0.4;
  }
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};



 