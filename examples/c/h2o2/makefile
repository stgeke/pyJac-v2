CXX = mpic++
CXXFLAGS= -fopenmp -O2 -ffast-math -march=native -mtune=native
DEPS =
LDFLAGS = -fopenmp 
OBJ = jacobian_main.o read_initial_conditions.o error_check.o timer.o jacobian_driver.o jacobian.o species_rates.o chem_utils.o

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
