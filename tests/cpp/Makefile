CXX ?= g++
PREFIX ?= $(CONDA_PREFIX)

CPPFLAGS := -I$(PREFIX)/include
CXXFLAGS := -O2 -g -std=c++14
LDFLAGS := -L$(PREFIX)/lib
LDLIBS := -larmadillo

ST_BIN := phyloacc_cpp_tests_st
GT_BIN := phyloacc_cpp_tests_gt
FAIL_ST_BIN := phyloacc_cpp_fail_st
FAIL_GT_BIN := phyloacc_cpp_fail_gt
FAIL_ST_PROFILE_BIN := phyloacc_cpp_fail_st_profile
FAIL_GT_PROFILE_BIN := phyloacc_cpp_fail_gt_profile

.PHONY: all clean

all: $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN)

$(ST_BIN): test_st_main.cpp ../../src/PhyloAcc-ST/newick.cpp ../../src/PhyloAcc-ST/profile.cpp ../../src/PhyloAcc-ST/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-ST $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(GT_BIN): test_gt_main.cpp ../../src/PhyloAcc-GT/newick2.cpp ../../src/PhyloAcc-GT/profile.cpp ../../src/PhyloAcc-GT/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-GT $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_BIN): fail_st_main.cpp ../../src/PhyloAcc-ST/newick.cpp ../../src/PhyloAcc-ST/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-ST $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_BIN): fail_gt_main.cpp ../../src/PhyloAcc-GT/newick2.cpp ../../src/PhyloAcc-GT/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-GT $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_PROFILE_BIN): fail_st_profile_main.cpp ../../src/PhyloAcc-ST/profile.cpp ../../src/PhyloAcc-ST/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-ST $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_PROFILE_BIN): fail_gt_profile_main.cpp ../../src/PhyloAcc-GT/profile.cpp ../../src/PhyloAcc-GT/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I../../src/PhyloAcc-GT $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN)
