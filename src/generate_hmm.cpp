/* 
 * File:                        generate_hmm.cpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        30/10/15
 * 1st operating system:        Linux [Ubuntu 3.13.0-37-generic]
 * 2nd operating system:        Mac OS X [El Capitan 10.11]       
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.6          
 */


#ifndef __GENERATE_HMM_CPP__
#define __GENERATE_HMM_CPP__

#include <iostream>
#include "../include/HMMGenerator.hpp"

int main (int argc, char* argv[])
{
   if (argc == 3) {

      // HMM generator instance that will compute probabilities without smoothing
      // HMMGenerator hmm_generator(argv[1],argv[2],false);

      // HMM generator instance that will perform standard add-one smoothing
      HMMGenerator hmm_generator(argv[1], argv[2]);

      // HMM generator instance that will perform smoothing with own specified summand
      // HMMGenerator hmm_generator(argv[1],argv[2],true,0.75);
   }

   else {
      std::cerr <<
              "USAGE: generate_hmm.cpp training_file output_file \n"
              "The training file has to contain one word and its corresponding pos tag "
              "per line, seperated by tabulators or spaces. Each beginning and end of a "
              "sentence in the data must be marked with an additional <BOS> <BOS> and "
              "respectively <EOS> <EOS> pair also seperated by tabulators. \n"
              "Some example files are included in /data. \n";
      exit(1);
   }
   return 0;
}

#endif