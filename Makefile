CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

COMMON_SRCS = fpi/fixed_point_iteration.cpp

MAIN_SRCS = main.cpp $(COMMON_SRCS)
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

all: solver

solver: $(MAIN_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN_OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(MAIN_OBJS) $(PERF_OBJS) solver Flux output.txt