#ifndef MORTON_CODER_HPP
#define MORTON_CODER_HPP
/** @file MortonCoder.hpp
 * @brief Define the MortonCoder class for Z-order-curve values, aka Morton
 *   codes.
 */

#include <cstdint>
#include <climits>
#include <cassert>

#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"


namespace detail {

/** Spreads the first 10-bits of a 32-bit number so that there are two 0s
 *   in between each bit.
 * @param[in] x 32-bit integer of the form 0b0000000000000000000000ABCDEFGHIJ
 * @returns     32-bit integer of the form 0b0000A00B00C00D00E00F00G00H00I00J,
 *   where each A,...,J is the corresponding bit from @a x
 */
inline uint32_t spread_bits(uint32_t x) {
  // HW4: YOUR CODE HERE
  x = (x | (x << 16)) & 0b0011000000000000000011111111;
  x = (x | (x <<  8)) & 0b0011000000001111000000001111;
  x = (x | (x <<  4)) & 0b0011000011000011000011000011;
  x = (x | (x <<  2)) & 0b1001001001001001001001001001;
  return x;
}

/** Compacts every third bit in a 32-bit integer to the lowest 10-bits.
 *   The inverse of spread_bits.
 * @param[in] x 32-bit integer of form 0bXXXXAXXBXXCXXDXXEXXFXXGXXHXXIXXJ
 * @returns     32-bit integer of form 0b0000000000000000000000ABCDEFGHIJ
 *   where each A,...,J is the corresponding bit from @a x
 */
inline uint32_t compact_bits(uint32_t x) {
  x &= 0b1001001001001001001001001001;
  x = (x | (x >> 2)) & 0b0011000011000011000011000011;
  x = (x | (x >> 4)) & 0b0011000000001111000000001111;
  x = (x | (x >> 8)) & 0b0011000000000000000011111111;
  x = (x | (x >> 2)) & 0b001111111111;
  return x;
}

/** Smears the bits in c into the low bits by steps of one
 *
 * Example: 00011100100 -> 000111111111
 */
inline uint32_t smear_low_1(uint32_t c) {
  c |= c >>  1;
  c |= c >>  2;
  c |= c >>  4;
  c |= c >>  8;
  c |= c >> 16;
  return c;
}

/** Smears the bits in c into the low bits by steps of three
 *
 * Example: 0000010000000000 -> 0000010010010010
 */
inline uint32_t smear_low_3(uint32_t c) {
  c |= c >>  3;
  c |= c >>  6;
  c |= c >> 12;
  c |= c >> 24;
  return c;
}

} // end namespace detail


/** @class MortonCoder
 * @brief Class representing Z-order-curve values, aka Morton codes.
 *
 * The Z-order curve is a space-filling curve: a one-dimensional curve that
 * fills a multi-dimensional space. Space-filling curves offer advantages for
 * representing points in 3D space. Points near each other in 3D space are
 * likely to have close Morton codes! So we can store points in a map with
 * Z-order value as key, then iterate over nearby Z-order values to fetch
 * points nearby in space.
 *
 * Unfortunately, it's impossible to reduce 3D space to 1D space perfectly:
 * there are some discontinuities in any space-filling curve mapping. But the
 * mapping is still an awesome tool, and some optimizations (see BIGMIN on the
 * Wikipedia page linked below) make them very effective.
 *
 * The MortonCoder class encapsulates a BoundingBox and can be used to translate
 * between spatial Points and Morton codes relative to that BoundingBox.
 *
 * A single Morton code corresponds to a rectangular volume within that
 * BoundingBox, called its <em>cell</em>. Each side of the BoundingBox is
 * divided into 2^L equal-sized cells, for a total of 8^L cells.
 *
 * Read more about the Z-order curve here:
 * http://en.wikipedia.org/wiki/Z-order_curve
 *
 * This class computes maps box numbers to point and visa-versa
 * with respect to a bounding box and the number of equal-volume boxes (8^L).
 * These mappings are performed in O(1) time.
 */
template <int L = 5>
class MortonCoder
{
 public:
  /** The type to use for the Morton codes -- allows 30-bit codes */
  using code_type = uint32_t;

  // Using a 32-bit unsigned int for the code_type
  // means we can only resolve 10 3D levels
  static_assert(L >= 1 && L <= 8*sizeof(code_type)/3,
                "L (LEVELS) must fit into code_type");

  /** The number of bits per dimension [octree subdivisions]. #cells = 8^L. */
  static constexpr int levels = L;
  /** The number of cells per side of the bounding box (2^L). */
  static constexpr code_type cells_per_side = code_type(1) << L;
  /** One more than the largest code (8^L). */
  static constexpr code_type end_code = code_type(1) << (3*L);

  /** Construct a MortonCoder with a bounding box. */
  MortonCoder(const Box3D& bb)
    : pmin_(bb.min()),
      cell_size_((bb.max() - bb.min()) / cells_per_side) {
    assert(!bb.empty());
  }

  /** Return the MortonCoder's bounding box. */
  Box3D bounding_box() const {
    return Box3D(pmin_, pmin_ + (cell_size_ * cells_per_side));
  }

  /** Return the bounding box of the cell with Morton code @a c.
   * @pre c < end_code
   */
  Box3D cell(code_type c) const {
    assert(c < end_code);
    Point p = deinterleave(c);
    p *= cell_size_;
    p += pmin_;
    return Box3D(p, p + cell_size_);
  }

  /** Return the Morton code of Point @a p.
   * @pre bounding_box().contains(@a p)
   * @post cell(result).contains(@a p)
   */
  code_type code(Point p) const {
    assert(bounding_box().contains(p));
    p -= pmin_;
    p /= cell_size_;
    return interleave(p);
  }

  // 0x09249249 = 0b001001001001001001001001001001
  static constexpr code_type coordinate_mask = 0x09249249;
  static constexpr code_type x_mask = coordinate_mask << 0;
  static constexpr code_type y_mask = coordinate_mask << 1;
  static constexpr code_type z_mask = coordinate_mask << 2;

  /** True if min <= idx <= max and idx is inside the box defined
   * by the Morton codes @a min and @a max
   */
  static bool is_in_box(code_type idx, code_type min, code_type max) {
    return (min & x_mask) <= (idx & x_mask) && (idx & x_mask) <= (max & x_mask)
        && (min & y_mask) <= (idx & y_mask) && (idx & y_mask) <= (max & y_mask)
        && (min & z_mask) <= (idx & z_mask) && (idx & z_mask) <= (max & z_mask);
  }

  /** Advance idx to the next box contained in the bounding box defined
   * by the Morton codes @a min and @a max
   *
   * @return result = @a min if @a idx <= @a min
   *                  @a idx if @a idx >= @a max
   *                  smallest i >= @a idx such that is_in_box(i,@a min,@a max)
   *                     and for all j in [@a idx,i), !is_in_box(j,@a min,@a max)
   * @post result > @a max || is_in_box(result,@a min,@a max)
   * @post result >= @a idx
   */
  static code_type advance_to_box(code_type idx, code_type min, code_type max) {
    if (idx >= max) return idx;

    // If outside the box in some coord, record the difference
    code_type delta = 0;
    if      ((idx & x_mask) > (max & x_mask))  delta |= (idx ^ max) & x_mask;
    else if ((idx & x_mask) < (min & x_mask))  delta |= (idx ^ min) & x_mask;
    if      ((idx & y_mask) > (max & y_mask))  delta |= (idx ^ max) & y_mask;
    else if ((idx & y_mask) < (min & y_mask))  delta |= (idx ^ min) & y_mask;
    if      ((idx & z_mask) > (max & z_mask))  delta |= (idx ^ max) & z_mask;
    else if ((idx & z_mask) < (min & z_mask))  delta |= (idx ^ min) & z_mask;

    // Delta is only zero if idx is in the box
    if (delta == 0) return idx;

    // Smear into a low bit mask, i.e. 0000111111111111
    delta = detail::smear_low_1(delta >> 1);
    // Chi masks high bits we cannot carry into
    const code_type chi = detail::smear_low_3(idx ^ max);
    // The first 0 in idx and chi that is higher than delta
    delta = ~(idx | ~chi | delta);
    delta = (delta & -delta) - 1;

    // Flip zero bit of idx and zero all lower bits
    idx = (idx | delta) + 1;

    // For each coordinate, if idx is low set to min
    if ((idx & x_mask) < (min & x_mask))  idx |= min & x_mask;
    if ((idx & y_mask) < (min & y_mask))  idx |= min & y_mask;
    if ((idx & z_mask) < (min & z_mask))  idx |= min & z_mask;

    return idx;
  }

 private:

  /** The minimum of the MortonCoder bounding box. */
  Point pmin_;
  /** The extent of a single cell. */
  Point cell_size_;

  /** Interleave the bits of n into x, y, and z.
   * @pre x = [... x_2 x_1 x_0]
   * @pre y = [... y_2 y_1 y_0]
   * @pre z = [... z_2 z_1 z_0]
   * @post n = [... z_1 y_1 x_1 z_0 y_0 x_0]
   */
  static inline code_type interleave(const Point& p) {
    return detail::spread_bits((code_type) p.x)
        | (detail::spread_bits((code_type) p.y) << 1)
        | (detail::spread_bits((code_type) p.z) << 2);
  }

  /** Deinterleave the bits from n into a Point.
   * @pre n = [... n_2 n_1 n_0]
   * @post result.x = [... n_6 n_3 n_0]
   * @post result.y = [... n_7 n_4 n_1]
   * @post result.z = [... n_8 n_5 n_2]
   */
  static inline Point deinterleave(code_type c) {
    return Point(detail::compact_bits(c),
                 detail::compact_bits(c >> 1),
                 detail::compact_bits(c >> 2));
  }
};

#endif
