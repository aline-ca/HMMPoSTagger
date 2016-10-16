
#################################################
#  Makefile (Linux/Mac OS) for HMM-PoS-Tagger.  #
#################################################

# 1. Type "make" in command line to build program.
# 2. To use g++ as compiler, simply replace it with clang++ as CPPCOMPILER.

CPPCOMPILER         = clang++
COMPILER_FLAGS      = -O3 -std=c++11
BOOST_INCLUDE	    = -I boost_1_59_0/

# Generate programs and documentation
all: build documentation

# Create executables
build: generate_hmm compute_forward_prob compute_viterbi
	
# Generate HMM file 
generate_hmm: src/generate_hmm.cpp include/HMMGenerator.hpp
	$(CPPCOMPILER) $(BOOST_INCLUDE) $(COMPILER_FLAGS) src/generate_hmm.cpp -o bin/generate_hmm

# Compute most likely tag sequence
compute_viterbi: src/compute_viterbi.cpp include/HMM.hpp
	$(CPPCOMPILER) $(BOOST_INCLUDE) $(COMPILER_FLAGS) src/compute_viterbi.cpp -o bin/compute_viterbi

# Compute forward probability
compute_forward_prob: src/compute_forward_prob.cpp include/HMM.hpp
	$(CPPCOMPILER) $(BOOST_INCLUDE) $(COMPILER_FLAGS) src/compute_forward_prob.cpp -o bin/compute_forward_prob

# Test target for demonstrating the forward algorithm (simply write "make test_forward" in command line)
test_forward: 
	bin/compute_forward_prob models/unicorn_summand_smoothed.hmm data/unicorn_test.txt
	
# Test target for demonstrating the viterbi algorithm (simply write "make test_viterbi" in command line)
test_viterbi: 
	bin/compute_viterbi models/unicorn_summand_smoothed.hmm data/unicorn_test.txt
	
# Generate documentation
documentation:
	doxygen doc/hmm_pos_tagger.doxygen
	mv html doc/

# Clean up (optional)
clean:
	rm bin/generate_hmm
	rm bin/compute_viterbi
	rm bin/compute_forward_prob
	rm -R doc/html