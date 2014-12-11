// xtree.cc

#include <vector>
#include <sstream>
#include "sym.h"
#include "xtree.h"
#include "utility.h"
#include <iostream>
#include <fstream>
#ifndef NOSERIALIZE
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif
template <typename tree_type>
void test() {
  tree_type* x1 = new tree_type("L1");
  tree_type* x2 = new tree_type("L2");
  tree_type* x3 = new tree_type("L3");
  tree_type* x4 = new tree_type("L4");
  tree_type* x5 = new tree_type("L5");

  x1->children.push_back(x2);
  x1->children.push_back(x3);
  x2->children.push_back(x4);
  x2->children.push_back(x5);

  std::cout << "*x1 = " << *x1 << std::endl;

  std::vector<symbol> ts;

  x1->terminals(ts);
  std::cout << "ts = " << ts << std::endl;

  delete x1;
}

//bool catcountyldtree_type<std::string>::compact_trees = false;
void test2() {
  typedef catcountyldtree_type tree_type;
  typedef catcountyldtree_type::yld_type yld_type;
  tree_type* x1 = new tree_type("L1",yld_type({ {symbol("a"), symbol("b"), symbol("c")}}));
  tree_type* x2 = new tree_type("L2",yld_type({ {symbol("a")} }));
  tree_type* x3 = new tree_type("L3");
  tree_type* x4 = new tree_type("L4");
  tree_type* x5 = new tree_type("L5");

  x1->children.push_back(x2);
  x1->children.push_back(x3);
  x2->children.push_back(x4);
  x2->children.push_back(x5);
  x1->count = 1;

  std::cout << "*x1 = " << *x1 << std::endl;

  const yld_type& ts  = x1->terminals();
  std::cout << "ts = " << ts << std::endl;

  tree_type* y1 = new tree_type("L1",yld_type({ {symbol("a"), symbol("b"), symbol("c")}}));
  tree_type* y2 = new tree_type("L2",yld_type({ {symbol("a")} }));
  tree_type* y3 = new tree_type("L3");
  tree_type* y4 = new tree_type("L4");
  tree_type* y5 = new tree_type("L5");

  y1->count = 1;
  y1->children.push_back(y3);//swap order of 2 and 3 wrt above
  y1->children.push_back(y2);
  y2->children.push_back(y5);
  assert(!(*y1 == *x1));
  y2->children.push_back(y4);
 
  std::cout << "*y1 = " << *y1 << std::endl;
  assert(*y1 == *x1);
  assert(*y2 == *x2);
#ifndef NOSERIALIZE
  std::ofstream ofs("arf");
  { boost::archive::text_oarchive oa(ofs);
    oa << y1;
  }
#endif
  y1->delete_tree();

#ifndef NOSERIALIZE
  std::ifstream ifs("arf");
  tree_type* read;
  { boost::archive::text_iarchive ia(ifs);
    ia >> read;
  }
  std::cout << "read= " << *read << std::endl;
  assert(*read == *x1);
  read->delete_tree();
#endif
  x1->delete_tree();
}

int main(int argc, char** argv) {
  test<cattree_type>();
  test<catcounttree_type>();
  test2();

}  // main()
