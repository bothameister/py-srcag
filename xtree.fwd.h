// Making this a template argument of catcountyldtree_type was too involved,
// so hard-set it here instead

#ifndef XTREE_FWD_H
#define XTREE_FWD_H
#include <vector>
#include <iostream>
#include <string>
#include "sym.h"
namespace {
typedef symbol S;
typedef std::vector<S> Ss;
typedef std::vector<Ss> Sss;
typedef Sss yld_type; 
}

std::ostream& operator<<(std::ostream& os, const ::yld_type& yld) {
  static const std::string bigsep =   "__"; //goes between ranges
  static const std::string smallsep = "_"; //goes between elements of a single range
  int i=0;
  for (const auto& a : yld) {
    int j=0;
    os << (i++ == 0 ? "" : bigsep);
    for (const auto& b : a)
      os << (j++ == 0 ? "" : smallsep) << b;
  }
  return os;
}

#endif /* end of include guard: XTREE_FWD_H */
