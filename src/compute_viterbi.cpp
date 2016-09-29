/* 
 * File:                        compute_viterbi.cpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        30/10/15
 * 1st operating system:        Linux [Ubuntu 3.13.0-37-generic]
 * 2nd operating system:        Mac OS X [El Capitan 10.11]       
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.6          
 */


#ifndef __COMPUTE_VITERBI_CPP__
#define __COMPUTE_VITERBI_CPP__

#include <iostream>
#include "../include/HMM.hpp"

// TODO 1: Fix debug switch/outputs.


int main (int argc, char* argv[])
{
  if (argc == 3) {

    /// HMM instance constructed from file in command line.
    /// Debug mode is set to false by default.
    HMM hmm(argv[1]);
    
    /// Set debugging parameter to true to enable debugging.
    // HMM hmm(argv[1], true);

    // Computes most likely sequence of hidden states for observation sequence 
    // read in from file in command line
    hmm.get_most_likely_seq(argv[2]);
    
  }

  else {
    std::cerr << 
      "USAGE: compute_viterbi hmm_file test_file\n\n"
      "hmm_file: File that contains the HMM. An example HMM file is included in /examples.\n"
      "test_file: File that contains test sentences for which the most likely tag sequence shall be computed.\n\n";
    exit(1);
  }
  return 0;
}

#endif 
