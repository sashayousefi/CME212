/**
 * @file mtl_test.cpp
 * Implimentation file for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// Define a IdentityMatrix that interfaces with MTL
class IdentityMatrix{
/** Helper function to perform multiplication . Allows for delayed
* evaluation of results .
* Assign::apply (a, b) resolves to an assignment operation such as
* a += b, a -= b, or a = b .
* @pre @a size (v) == size (w) */
    public:
    int s;
    //constructor for identitymatrix
    IdentityMatrix(int _s) : s(_s * _s) {}
    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn& v, VectorOut& w, Assign) const{
        Assign::apply(w, v);
    }
    /* * Matvec forwards to MTL ’s lazy mat_cvec_multiplier operator */
    template <typename Vector>
    mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
    operator*(const Vector& v) const {
        return {*this, v};
    }
    private:
    // Empty!
};
/* * The number of elements in the matrix . */
inline std::size_t size(const IdentityMatrix& A){
    return A.s * A.s;
}

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const IdentityMatrix& A){
    return A.s;
}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const IdentityMatrix& A){
    return A.s;
}
