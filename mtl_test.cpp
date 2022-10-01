#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include "IdentityMatrix.hpp"

/** Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl {
namespace ashape {

/** Define IdentityMatrix to be a non - scalar type . */
template <>
struct ashape_aux<IdentityMatrix> {
    typedef nonscal type;
    };
} // end namespace ashape

/** IdentityMatrix implements the Collection concept
* with value_type and size_type */
template <>
struct Collection<IdentityMatrix> {
    typedef double value_type;
    typedef unsigned size_type;
    };
} // end namespace mtl


int main()
{
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
  const int size = 10, N = size * size; 

  typedef IdentityMatrix matrix_type;
  // Set up a matrix
  matrix_type I(size);

  // Create an preconditioner
  itl::pc::identity<matrix_type> L(I);

  // Set b such that x == 1 is solution; start with x == 0
  mtl::dense_vector<double> x(N, 1.0), b(N);
  b = I*x; x=0;
  assert(b != x);
  
  // Termination criterion: r < 1e-10 * b or N iterations
  itl::noisy_iteration<double> iter(b, 500, 1.e-6);

  // Solve Ax == b with left preconditioner P
  itl::cg(I, x, b, L, iter);

  return 0;
}
