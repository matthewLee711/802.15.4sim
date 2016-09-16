###################################################################
#  Makefile for example matrix class
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  2 July 2015
###################################################################

# compiler & flags
CXX = g++
CXXFLAGS = -O -std=c++11
#CXXFLAGS = -O0 -g -std=c++11

# makefile targets
all : matrix_test.exe intro.exe

matrix_test.exe : matrix_test.cpp matrix.o GramSchmidt.o
	$(CXX) $(CXXFLAGS) $^ -o $@

intro.exe : intro.cpp matrix.o
	$(CXX) $(CXXFLAGS) $^ -o $@

matrix.o : matrix.cpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

GramSchmidt.o : GramSchmidt.cpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean :
	\rm -f *.o *.out a_data *.txt

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
