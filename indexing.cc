#include <cstdlib>
#include <iostream>
#include <cassert>
#include <set>
#include <functional>
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include "nested_list.h"
#include "consts.h"

using namespace std;

typedef mcfg::U U;
typedef mcfg::R R;
typedef mcfg::rvec_t rvec_t;

typedef std::initializer_list<U> rvec_lt;
typedef std::unordered_map<R,rvec_t> setR;
int debug=10;
bool verbose = true;
//#define get_MAX_INPUT_LEN 20 //uncomment to overwrite the one in consts.h

void verify(std::initializer_list<U> rv, setR& tracker, bool* flag=0)
{
  auto id = mcfg::index_bitvec<U,R,get_MAX_INPUT_LEN>(rv);
  rvec_t proper_rv;
  list2rvec(rv, proper_rv);
  if (verbose)
    cout << proper_rv << " => " << id << std::endl;
  if (tracker.find(id) != tracker.end()) {
    cout << proper_rv << " clashes with " << tracker[id];
    cout << ". Key=" << id << endl;
    if (flag)
      *flag=false;
    else
      assert(false);
  }
  else
    tracker.insert(make_pair(id, proper_rv));
}

//The structure of this method should mimick that of the actual inside algorithm in py-cky.h
//for this testing to be of use.
//Returns true if no clashes
bool verify(U n, int max_rangevector_dim) {
  setR track;
  bool flag=true;
  bool* fp = &flag;

  for (U i=0;i<n; ++i)
    verify({i,i+1}, track, fp);

  for (U gap = 2; gap <= n; ++gap) { // non-terminals. e.g. rule shape A -> B C
    for (U left = 0; left + gap <= n; ++left) {
      U right = left + gap;

      for (U c = left+1; c < right && max_rangevector_dim >= 2 ; ++c) { // c and d are indexing 
        for (U d = c; d < right; ++d) {                                 // for parent_ranges_2 

          for (U e = d+1; e < right && max_rangevector_dim >= 3; ++e) { // e and f are indexing
            for (U f = e; f < right; ++f) {                             // for parent_ranges_3
             
              for (U g = f+1; g < right && max_rangevector_dim >= 4; ++g) { // g and h are indexing
                for (U h = g; h < right; ++h) {                             // for parent_ranges_4

                  for (U i = h+1; i < right && max_rangevector_dim >= 5; ++i) { // i and j are indexing
                    for (U j = i; j < right; ++j) {                             // for parent_ranges_5
                      //Deal with cases where parent rangevec has 5 ranges
                      verify({left,c,d,e,f,g,h,i,j,right}, track, fp);
                    }
                  }
                  //Deal with cases where parent rangevec has 4 ranges
                  verify({left,c,d,e,f,g,h,right}, track, fp);
                }
              }
              //Deal with cases where parent rangevec has 3 ranges
              verify({left,c,d,e,f,right}, track, fp);
            }
          }

          //Deal with cases where parent rangevec has 2 ranges
          verify({left,c,d,right}, track, fp);
        }
      }
      //Deal with cases where parent rangevec has 1 range
      verify({left,right}, track, fp);
    }
  }
  return flag;
}



int main(int argc, const char *argv[])
{ 
  verbose= (argc>1);
  double bits = log(numeric_limits<R>::max())/log(2);
  cout << "Test for chart indexing and bit-space\n";
  cout << "=====================================\n";

  cout << "Chart is indexed with type R, which is " << bits << "-bit on this machine." << endl;
  cout << "At maximum input length ...., the implementation can handle up to .... ranges per rangevector.\n";
  
  for (int i=5; i<= get_MAX_INPUT_LEN; i+=5) {
    int blocks_avail =floor(bits / (ceil(log(i+1)/log(2)))); 
    cout << "\t\t" << i << "\t\t" << blocks_avail/2 << endl;
    //cout << "   at max input len " << i << " impl can handle up to " << blocks_avail/2 << " ranges in a rangevec\n";
  }
  cout << "\nTesting for key collisions...\n";
  for (int i=get_MAX_INPUT_LEN; i>0; i-=3)
  {
    cout << "n=" << i << "-----------\n";
    if (verify(i, get_MAX_RANGEVECTOR_DIM_IMPL))
      cout << "\tOK\n";
    else {
      cout << "\tBad.\n";
      //break;
    }
  }



  return 0;
}
