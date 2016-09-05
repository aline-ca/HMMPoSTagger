
#################################################
#  Makefile (Linux/Mac OS) for HMM-PoS-Tagger.  #
#################################################

# 1. Type "make" in command line to build program.
# 2. To use g++ as compiler, simply replace it with clang++ as CPPCOMPILER.
# 3. Please change the path in BOOST_INCLUDE to your system's boost library path. 

CPPCOMPILER         = clang++
COMPILER_FLAGS      = -O3 -std=c++11
BOOST_INCLUDE	    = -I /home/aline/boost_1_59_0/

# Generate programs and documentation
all: build

# Create executables
build: generate_hmm main

# Generate HMM file 
generate_hmm: src/generate_hmm.cpp include/HMMGenerator.hpp
	$(CPPCOMPILER) $(BOOST_INCLUDE) $(COMPILER_FLAGS) src/generate_hmm.cpp -o bin/generate_hmm

# main file
main: src/main.cpp include/HMM.hpp
	$(CPPCOMPILER) $(BOOST_INCLUDE) $(COMPILER_FLAGS) src/main.cpp -o bin/main


# Clean up (optional)
clean:
	rm bin/generate_hmm
	rm bin/main
