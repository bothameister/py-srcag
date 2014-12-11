# Modified: Jan Botha (April 2013)
# Original author: Mark Johnson (27th August 2009)
#
# Configuration:
BOOST_LIBDIR=/home/jan/loro/only_boost155/lib
BOOST_INCLUDEDIR=/home/jan/loro/only_boost155/include

# Usage:
#   production:
#      $ make clean py-srcag yieldextractor
#   debugging:
#      $ make clean run-toy NDEBUG=
#
#
# To build without needing boost, add NOSERIALIZE=-DNOSERIALIZE to the command lines above.

TARGETS := py-srcag yieldextractor indexing

top: $(TARGETS)

.PHONY: run-toy
run-toy: py-srcag testengger.mlt testeng.yld
	$(EXEC) ./py-srcag -F trace -D -R -1 -d 500 -a 0 -b 1 -w 1 -e 1 -f 1 -g 10 -h 0.01 -n 1000 -E -A testeng.prs -N 10 -G testeng.wlt testengger.mlt  < testeng.yld

testengger.mlt:
	./preprocess-grammar.py < testengger.lt > testengger.mlt

.PHONY: run-toyx
run-toyx: py-srcag testengger.mltx testeng.yld
	$(EXEC) ./py-srcag -F trace -D -R -1 -d 500 -a 0 -b 1 -w 1 -e 1 -f 1 -g 10 -h 0.01 -n 1000 -E -A testeng.prs -N 10 -G testeng.wlt testengger.mltx  < testeng.yld
############################################################
#                                                          #
#                    Program build                         #
#                                                          #
############################################################

# MODES: PRODFAST = fast & no assertions; PROD = fast with assertions; DEBUG = no optimization
MODE=PROD

ifeq ($(MODE),PROD)
  #
  # production
  #
  CC = $(CXX)
  NDEBUG=
  CFLAGS = -MMD -O6 -Wall -ffast-math -finline-functions -fomit-frame-pointer -fstrict-aliasing $(GCCFLAGS)
  LDFLAGS = -Wall -O6 $(GCCLDFLAGS)
  EXEC = time
else ifeq ($(MODE),PRODFAST)
  CC = $(CXX)
  NDEBUG=-DNDEBUG
  CFLAGS = -MMD -O6 -Wall -ffast-math -finline-functions -fomit-frame-pointer -fstrict-aliasing $(GCCFLAGS)
  LDFLAGS = -Wall -O6 $(GCCLDFLAGS)
else
  #
  # debugging
  #
  NDEBUG=
  CFLAGS = -g -O0 -MMD -Wall -ffast-math -fstrict-aliasing $(GCCFLAGS)
  LDFLAGS = -g -Wall $(GCCLDFLAGS)
  EXEC = valgrind --tool=memcheck --leak-check=full
endif

NOSERIALIZE=
EMPTY=
ifeq ($(NOSERIALIZE), $(EMPTY))
  BOOST_SERIALIZE_LDFLAGS = -L$(BOOST_LIBDIR) -lboost_serialization 
  BOOST_INCLUDE = -I$(BOOST_INCLUDEDIR)
else
  NOSERIALIZE = -DNOSERIALIZE
  BOOST_SERIALIZE_LDFLAGS =
  BOOST_INCLUDE =
endif
#
# profiling
#
# CFLAGS = -g -pg -O6 -MMD -Wall -ffast-math -fno-default-inline -fno-inline $(GCCFLAGS)
# CFLAGS = -g -pg -O -MMD -Wall -ffast-math $(GCCFLAGS)
# LDFLAGS = -g -pg

LDFLAGS += $(BOOST_SERIALIZE_LDFLAGS)
CXXFLAGS = $(CFLAGS) $(NOSERIALIZE) $(NDEBUG) -std=c++0x -I. $(BOOST_INCLUDE)
SOURCES = gammadist.c mt19937ar.c py-srcag.cc xtree.cc sym.cc yieldextractor.cc
OBJECTS = $(patsubst %.l,%.o,$(patsubst %.c,%.o,$(SOURCES:%.cc=%.o)))

py-srcag: gammadist.o py-srcag.o mt19937ar.o sym.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

yieldextractor: gammadist.o yieldextractor.o mt19937ar.o sym.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

indexing: indexing.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

gammadist.o: gammadist.c
	gcc -c $(CFLAGS) -std=c99 $< -o $@

mt19937ar.o: mt19937ar.c
	gcc -c $(CFLAGS) $< -o $@

test: consistency.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< && ./$@

treetest: xtree.o sym.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) && ./$@

xtree.o: xtree.cc xtree.h xtree.fwd.h
	$(CXX) $(CXXFLAGS) -o $@ -c xtree.cc


.PHONY: 
clean: 
	rm -fr *.o *.d *~ core test

.PHONY: real-clean
real-clean: clean 
	rm -fr $(TARGETS)

# this command tells GNU make to look for dependencies in *.d files
-include $(patsubst %.l,%.d,$(patsubst %.c,%.d,$(SOURCES:%.cc=%.d)))
