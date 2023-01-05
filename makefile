CXX := g++
CXXFLAGS := -O3 -Wall -Wextra -Wconversion -std=c++11

TARGET := main
OBJS := main.o COO2CSR.o CGSolver.o matvecops.o heat.o sparse.o
INCS := COO2CSR.hpp CGSolver.hpp matvecops.hpp heat.hpp sparse.hpp

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS)

%.o: %.cpp $(INCS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

.PHONY: clean
clean: 
	$(RM) $(OBJS) $(TARGET)




