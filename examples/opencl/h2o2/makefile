CXX = g++
CXXFLAGS=-x c++
DEPS =
LDFLAGS = -L/soft/compilers/cuda/cuda-10.0.130/lib64 -lOpenCL 
OBJ = jacobian_main.o read_initial_conditions.o error_check.o timer.o jacobian_compiler.o

%.o: %.ocl $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all clean

all: genic jac

jac: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

genic: genic.cpp
	$(CXX) -o $@ genic.cpp

clean:
	rm -f *.o jac genic jacobian.bin 
