#ifndef MCFG_H_4XYWFJVK
#define MCFG_H_4XYWFJVK

#include <iostream>
#include <vector>
#include <sstream>
#include <limits>
#include <cassert>
#include <boost/functional/hash.hpp>
#include "utility.h"
extern int debug;

#define sayif(msg) if (debug > 20500) std::cout << msg << std::endl
/////// Compile with -std=c++0x
namespace mcfg {

  //! Creates a hash of the values, which are interpreted as a sequence of ranges
template<typename U>
U index(const std::initializer_list<U>& values) {
  assert(values.size() % 2 == 0);
  return boost::hash_range(begin(values), end(values));
}

//Constructs a hash of values by reserving a block for each one and pushing its bits left
//U: type of elements in 'values'. R: return type (should be big type). L: max input length to support
template<typename U, typename R, unsigned L>
R index_bitvec(const U* values, U size) {
  static const float log2 = log(2);
  static const short bits_per_block = ceil(log(L+1)/log2); //num bits needed to store L. e.g. L='15' -> 4; '16' -> 5
  static const R preshiftmax = std::numeric_limits<R>::max() >> bits_per_block;

  assert(size % 2 == 0);
  R result=0;
  for (const U* p = values; p < values+size; ++p) {
    const auto& v = *p;
    assert(v <= L); //v has to fit into single block of bits
    if (result >= preshiftmax) {
      std::cerr << "fail! result >= preshiftmax: " << result << " >= " << preshiftmax << ", L=" << L << "bpb=" << bits_per_block << ", v=" << v << std::endl;
    }
    assert(result < preshiftmax); //result needs space for another block

    result <<= bits_per_block;
    result |= v;
  }
  return result;
}

//Alternative interface to the main index_bitvec, this one takes an initializer_list
template<typename U, typename R, unsigned L>
inline R index_bitvec(const std::initializer_list<U>& values) {
  return index_bitvec<U,R,L>(values.begin(), values.size());
}

// Convenience function to access a particular item in an std::initializer_list
template<typename T, typename I> 
inline const T&  at(const std::initializer_list<T>& list, I i)
{
  assert(list.size() > 0); //otherwise *begin(list) will segfault
  assert(i >= 0 && i < (I)list.size());
  return *(begin(list) + i);
}

typedef short Pos;
struct range_t {
  range_t(Pos left, Pos right) : left(left), right(right) { assert(sane()); }

  // Special marker values to denote 'empty' range; useful when having to allocate a range vector
  range_t() : left(std::numeric_limits<Pos>::max()), right(std::numeric_limits<Pos>::min()) {}

  //designates the range
  Pos left, right;
  
  static range_t empty;
  bool is_empty() const {
    return *this == empty;
  }

  //Comparators skip actual comparison when this/either object is 'empty'
  bool sane() const {
    return is_empty() || left < right;
  }

  friend bool operator==(const range_t& o1, const range_t& o2) {
    return o1.left == o2.left && o1.right == o2.right;
  }
  friend bool operator!=(const range_t& o1, const range_t& o2) {
    return !(o1 == o2); 
  }
  friend bool adjacent(const range_t& o1, const range_t& o2) {
    return either_empty(o1,o2) || 
               (o1.left < o2.left && o1.right == o2.left)   // o1o2
            || (o2.left < o1.left && o2.right == o1.left);  // o2o1
  }
  // o1 is left of o2 when o1 ends before at start of o2
  // i.e. two distinct ranges may be adjacent
  friend bool before(const range_t& o1, const range_t& o2) {
    return either_empty(o1,o2) || 
            o1.right <= o2.left;
  }
 
  static bool either_empty(const range_t& o1, const range_t& o2) {
    return  o1.is_empty() || o2.is_empty();
  }

  friend std::ostream& operator<<(std::ostream& os, const range_t& obj) {
    return os << std::make_pair(obj.left, obj.right);
  }
};
range_t range_t::empty;

typedef std::vector<range_t> rvec_t;
typedef std::vector<rvec_t> rvecs_t;

///// Types for the variable-to-range mapping
typedef int Idx;
typedef std::pair<Idx,Idx> IdxPair;
typedef std::vector<IdxPair> IdxPairs;
typedef std::vector<IdxPairs> Mapping;

typedef unsigned short int ushint;
typedef std::pair<ushint, ushint> intpair;
// returns the maximal dimension of the predicates involved in the mapping m
// pair.first is the max lhs dim; pair.second the max rhs dim.
// e.g. A(x,y) -> B(x) Y(y) gives 2,1
//      A(axyz) -> B(a,x,z) C(y) gives 1,3
intpair get_max_pred_dim(const Mapping& m) {
  ushint lhs_arg_max = 0;
  ushint rhs_arg_max = 0;
  for (const auto& rhs_vars : m) {
    if ((ushint)rhs_vars.size() > rhs_arg_max) rhs_arg_max = rhs_vars.size();
    for (const auto& var : rhs_vars)
      if ((ushint)var.first+1 > lhs_arg_max) lhs_arg_max = var.first+1; //+1 because var.first is a zero-based index
  }

  return std::make_pair(lhs_arg_max, rhs_arg_max);
}
// Figure out the dimensionality of the result vector subject to the given mapping.
// Instantiates empty result vector
rvecs_t instantiate_result_holder(const Mapping& m) {
  int lhs_arg_max = 0;
  for (const auto& rhs_vars : m)
    for (const auto& var : rhs_vars)
      if (var.first > lhs_arg_max) lhs_arg_max = var.first;

  rvecs_t lhs_vars;
  for (int lhs_arg_i = 0; lhs_arg_i <= lhs_arg_max; ++lhs_arg_i) {
    int this_arg_max = 0;
    for (const auto& rhs_vars : m)
      for (const auto& var : rhs_vars)
        if (var.first==lhs_arg_i && var.second > this_arg_max) this_arg_max = var.second;
    lhs_vars.push_back(rvec_t(this_arg_max+1, range_t()));
  }
  return lhs_vars;
}

// Determines whether the proposed compose action is allowed.
// The 'action' is to integrate the range vector passive into holder according to
// the mapping m. If successful, the caller can obtain the result in 'out'.
// Example of invalid cases:
//   - holder already has a range in a certain target position
//   - consecutive range vectors that need to form a single range later aren't consecutive
//   - ordering/overlap violations
bool compose(const rvecs_t& holder, const IdxPairs& m, const rvec_t& passive, rvecs_t& out) {
  assert(passive.size() > 0);
  assert(passive.size() == m.size());
  //copy so that we can modify out even if not sure if all actions will be allowed
  out.reserve(holder.size());
  std::copy(holder.begin(), holder.end(), std::back_inserter(out));
  for (unsigned i = 0; i < passive.size(); i++) {
    auto& range = passive[i];
    auto idx_pair = m[i]; //describes where in 'out' 'range' has to go
    
    auto& arg = out[idx_pair.first];
    auto& target_range = arg[idx_pair.second];
    Idx last_arg_idx = (Idx)arg.size() - 1; 
    //Format of these conditional statements are
    //  if <situation> && /**/ !<statement that should hold under this situation>
    //The ordering among conditions that overlap logically should be that the stricter comes first
    if (!target_range.is_empty()) {
      sayif("\tNO-occupied:" << target_range << ", " << range_t::empty); return false; }
    else 
    { 
      //if the lhs arg has a range before & after target, both sides have to be adjacent
      if (idx_pair.second > 0 && idx_pair.second < last_arg_idx 
                  && /**/ !(adjacent( arg[idx_pair.second-1], range)
                             && adjacent(range, arg[idx_pair.second+1]) ))
        {sayif("\tNO-non-adj"); return false; }

      //if the lhs arg has a range before target, they have to be adjacent/concatenate
      else if (idx_pair.second > 0 && /**/ !adjacent( arg[idx_pair.second-1], range))
        {sayif("\tNO-left-non-adj"); return false; }

      //if the lhs arg has a range after target, they have to be adjacent/concatenate
      else if(idx_pair.second < last_arg_idx  && /**/ !adjacent(range, arg[idx_pair.second+1] ))
        {sayif("\tNO-right-non-adj"); return false; }

      // if targeting the first range in an lhs arg, check consistency with end of preceding arg
      else if (idx_pair.second == 0 && idx_pair.first > 0 && /**/ !before(out[idx_pair.first-1].back(), range))
        {sayif("\tNO-prev-arg-ends-after-this"); return false; }

      // if targeting the last range in an lhs arg, check consistency with start of next arg
      else if (idx_pair.second == last_arg_idx && idx_pair.first < (Idx)out.size()-1 && /**/ !before(range, out[idx_pair.first+1].front())) 
        { sayif("\tNO-next-arg-starts-before-this"); return false; }

      else
      {
        target_range = range;
        sayif("\tOK, yields " << out);
      }
    }
  }
  return true;
}

template<typename U>
void list2rvec(const std::initializer_list<U>& rvec_values, rvec_t& rvec) {
  assert(rvec_values.size() % 2 == 0);
  for (auto x=begin(rvec_values); x<end(rvec_values);x+=2)
    rvec.push_back(range_t(*x, *(x+1)));
}

template<typename U>
bool result_dims_match(const rvecs_t& found, const std::initializer_list<U>& exp_) {
    rvec_t exp; list2rvec(exp_, exp);
    return exp.size() == found.size();
}

template<typename U>
bool compose_activepassive( const std::initializer_list<U>& active_rvec_values, const std::initializer_list<U>& passive_rvec_values, const Mapping& m,  const std::initializer_list<U>& expected_result)
{
  assert(m.size() == 2); //binarised grammar assumption
  
  if (active_rvec_values.size() / 2 != m[0].size() //fail early if dimensionality incompat.
      || passive_rvec_values.size() / 2 != m[1].size()) {
    sayif("\tNO, mismatch between rvec-dim and m-dim");
    return false;
  }

  rvecs_t o = instantiate_result_holder(m);
  if (!result_dims_match(o, expected_result)) { //fail early if dimensionality incompat.
    sayif("\tNO, mismatch between result-dim and expected-dim");
    return false;
  }
  rvec_t active_rvec; list2rvec(active_rvec_values, active_rvec);
  rvec_t passive_rvec; list2rvec(passive_rvec_values, passive_rvec);
  rvecs_t o1, o2;
  if (!compose(o, m[0], active_rvec, o1)) {
    sayif("\tNO, couldn't even place the first rvec into holder");
    return false;
  }
  return compose(o1,m[1], passive_rvec, o2); //o2 is potentially useful to the caller...
}

template<typename U>
bool compose_unary(const std::initializer_list<U>& child_rvec_values, const Mapping& m,  const std::initializer_list<U>& expected_result)
{
  assert(m.size() == 1); //only treat unary rules here
  const IdxPairs& m0 = m[0];
  if (child_rvec_values.size() / 2 != m0.size()) {
    sayif("\tNO, mismatch between rvec-dim and m-dim");
    return false;
  }

  rvecs_t o = instantiate_result_holder(m);
  if (!result_dims_match(o, expected_result)) { //fail early if dimensionality incompat.
    sayif("\tNO, mismatch between result-dim and expected-dim");
    return false;
  }
  rvec_t child_rvec; list2rvec(child_rvec_values, child_rvec);
  rvecs_t o1;
  bool possible = compose(o, m0, child_rvec, o1);
  if (!possible) return false;
  // Check for the following
  // expected: ( (0 1)        (1 2)           (2 4) )
  // o1      : (((0 1))   ((1 2) (2 3))         ((3 4))) << bad
  // o1'     : (((0 1))      ((1 2))      ((2 3) (3 4))) << OK
  // Through concatenation, it's possible to get 
  // an internally consistent composition o1
  // whose top-level dim matches expected_result, but
  // violates its range boundaries
  sayif("\t=> Possible, checking extents against expected result..."); 
  std::size_t i=0;
  for (auto x=begin(expected_result); x<end(expected_result);x+=2,++i) {
    Pos exp_left=*x, exp_right=*(x+1);
    if (exp_left!= o1[i].front().left || exp_right != o1[i].back().right)
      return false;
  }
  return true;
}
} /* namespace mcfg */

#endif /* end of include guard: MCFG_H_4XYWFJVK */
