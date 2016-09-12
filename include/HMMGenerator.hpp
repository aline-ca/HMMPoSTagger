/* 
 * File:                        HMMGenerator.hpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        30/10/15
 * 1st operating system:        Linux [Ubuntu 3.13.0-37-generic]
 * 2nd operating system:        Mac OS X [El Capitan 10.11]       
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.6          
 */

/* 
  TODO: 
    - Implement better smoothing algorithm (Good Turing, Witten-Bell or Absolute Discounting)
*/

#ifndef __HMM_GENERATOR_HPP__
#define __HMM_GENERATOR_HPP__

#include <boost/tokenizer.hpp>  // boost tokenizer
#include <iostream>             // std::cout
#include <cstdlib>              // exit()
#include <fstream>              // std::ifstream
#include <cassert>              // assert()
#include <algorithm>            // std::sort
#include <string>               // std::string
#include <tuple>                // std::tuple (requires C++11)
#include <vector>               // std::vector
#include <map>                  // std::map
#include <set>                  // std::set

using namespace boost;

/** 
  @brief The HMMGenerator class is used to read in a tsv file that contains words and their corresponding PoS-tags.
         It creates a HMM from it and writes the properties of the trained HMM into a new tsv file.
         The HMM class can then read in the HMM and do forward/viterbi calculations.   
  @note The implemented HMM is case sensitive, e.g. "The" and "the" are seen as two different words.
  @note Futhermore, all outgoing transitions for a state add up to one (roughly) EXCEPT for the transitions from EOS 
        (since it is the last word in the file).  
*/
class HMMGenerator
{
 public: // Types 
  bool                                        smoothing;          ///< Smoothing configuration. Set to true by default
  double                                      smoothing_summand;  ///< For smoothing computation. Set to one by default (Add-one)

 private: // Typedefs
  typedef std::vector<std::string>            StringVector;       ///< String vector for reading in tags and words
  typedef std::tuple<std::string,std::string> BigramTuple;        ///< String tuple as bigram representation
  typedef std::set<std::string>               StringSet;          ///< String set for all observed words.
  typedef std::map<std::string,int>           StringIntMap;       ///< String-to-int map for counting word occurrences
  typedef std::map<BigramTuple,int>           TupleIntMap;        ///< Tuple-to-int map for counting bigrams and observations
  typedef std::map<std::string,double>        StringDoubleMap;    ///< String-to-double map for saving initial probabilities
  typedef std::map<BigramTuple,double>        TupleDoubleMap;     ///< Tuple-to-double map for saving MLE and word likelihood
  typedef std::map<std::string,unsigned>      StringUnsignedMap;  ///< String-to-unsigned map for representing state as number
  typedef tokenizer< char_separator<char> >   Tokenizer;          ///< Boost tokenizer used for tsv file parsing

 private: // Types
  unsigned          token_counter;            ///< Number of all tokens in the corpus
  StringIntMap      type_counter_map;         ///< Maps a pos tag to the number of its occurrences
  StringVector      pos_seq;                  ///< Vector of all PoS tags in the order they occurred in the corpus
  StringVector      word_seq;                 ///< Vector of all words in the order they occurred in the corpus
  StringSet         observed_words_set;       ///< Set that saves all observed words
  TupleIntMap       bigram_counter_map;       ///< Maps a bigram to the number of its occurrences 
  TupleIntMap       observation_counter_map;  ///< Maps a tag-word-tuple to the number of its occurrences
  TupleDoubleMap    ml_map;                   ///< Saves maximum likelihood 
  TupleDoubleMap    word_likelihood_map;      ///< Saves word likelihood
  StringUnsignedMap tag_state_map;            ///< Maps a pos tag to its state representation
  StringDoubleMap   init_probs;               ///< Maps a pos tag to the probability to start a sentence with it

 public: // Functions

  /**
    @brief Constructor that reads in a tsv file containing words and their corresponding pos-tags, 
           does all necessary computations to create a HMM from it and outputs the HMM in tsv format.
    @param infile Pointer to file that will be read in
    @param outfile Pointer to file in which the output will be written
    @param summand Number that will be used as summand for smoothing computations. Is set to one by default.
    @param sm Bool that specifies whether smoothing should be enabled. Is set to true by default. 
  */
  HMMGenerator(const char* infile, const char* outfile, bool sm=true, double summand=1) {
    smoothing = sm;               // Set smoothing variable to specified configuration
    smoothing_summand = summand;  // Set specified number as smoothing summand
    token_counter = 0;            // Initialize token counter
    read_in_file(infile);         // Read in input file
    get_bigrams();                // Get all bigrams that occurred in the file
    remove_duplicates(pos_seq);   // Remove duplicates from pos_seq for further use
    compute_mle();                // Compute maximum likelihood estimation
    compute_word_likelihood();    // Compute word likelihood
    initialize_init_probs();      // Initialize initial probabilities map
    get_init_probs();             // Collect initial probabilities 
    build_tag_state_map();        // Create map that maps a pos tag to its state index
    generate_tsv_file(outfile);   // Output tsv file that contains the HMM
  }   


 private: // Functions
  /**
    @brief Reads in a file, tokenizes it and stores the information. Saves the
           sequence of occurred words, the sequence of occurred pos-tags, the 
           number of counts, the number of types and the number of observations. 
    @param filename Pointer to file that will be read from
  */
  void read_in_file(const char* filename) {
      
    if (filename != 0) {                        // If pointer is no null pointer 
      std::string line;                         // Current line in file
      StringVector vec;                         // Vector that saves the current line
      std::ifstream in(filename);               // Ifstream object constructed from file  
      
      if (in.is_open()) {                       // If it was possible to open the file
          
        while (getline(in,line)) {              // Read in each line seperately
          char_separator<char> sep("\t ");      // Define tab and space field separators    
          // Tokenizer constructed from current line and defined separator: 
          Tokenizer tok(line, sep);
          vec.assign(tok.begin(),tok.end());    // Assign current line to vector
          // Ignore all lines that do not have exactly two fields (word and pos tag):
          if (vec.size() != 2) continue;
          
          const std::string& word = vec[0];     // Current word
          const std::string& pos_tag = vec[1];  // Current pos tag
          
          word_seq.push_back(word);             // Append word to word sequence
          pos_seq.push_back(pos_tag);           // Append pos tag to pos sequence
          
          // Increase type and token counter by one:
          ++type_counter_map[pos_tag];
          ++token_counter;
          
          // Add observed word to the set:          
         if (observed_words_set.find(word) == observed_words_set.end()) {
             observed_words_set.insert(word); 
         }
          
          // Pointer to tuple that represents reading pos_tag and emitting word at the same time: 
          auto observation = std::make_tuple(pos_tag, word);  
          // Increment map counter at that position:
          ++observation_counter_map[observation];
        } 
      }
      else {
        std::cerr << "Unable to open '" << filename << "'\n";
        exit(1);
      }
    }
    else {
      std::cerr << "Memory error: Pointer to null. Could not read in files.\n";
      exit(1);
    }
  }


  /**
    @brief Generates a tsv file that contains all properties of the HMM 
           and can be further processed by the HMM class.
    @param filename Pointer to file that will be written in
  */
  void generate_tsv_file(const char* filename) const 
  {
      
    if (filename != 0) {              // If pointer is no null pointer
      std::ofstream out(filename);    // Ofstream object created from file
      
      if (out.is_open()) {            // If it was possible to open the file
          
        // Write metadata (num of states + num of observations)
        out << "#\t" << tag_state_map.size() << "\t" << observed_words_set.size() << "\n";
                
        // Write all transitions into file by iterating over bigram-probability map:
        out << "\nTransitions:\n";
        
        for (auto e = ml_map.begin(); e != ml_map.end(); ++e) {
          std::string tag1;           // First tag of bigram
          std::string tag2;           // Second tag of bigram
          std::tie (tag1, tag2) = e->first;
          out << get_state(tag1) << "\t" << get_state(tag2) << "\t" <<  e->second << "\n";
        } 
        out << "\n~~~~~~~~~~~~~~~~~~~~~~~~\n";

        // Write all state-symbol mappings into file by iterating over tag state map:
        out << "\nState Symbol Map:\n";
        
        for (auto e = tag_state_map.begin(); e != tag_state_map.end(); ++e) {
            out << e->second << "\t" << e->first << "\n"; 
        }
        out << "\n~~~~~~~~~~~~~~~~~~~~~~~~\n";

        // Write all observations into file by iterating over word likelihood map:
        out << "\nObservations:\n";
        
        for (auto e = word_likelihood_map.begin(); e != word_likelihood_map.end(); ++e) {
          std::string tag;        // Current tag
          std::string word;       // Current word
          std::tie (tag, word) = (e->first);
          out << tag << "\t" << word << "\t" << e->second << "\n";
        }
        out << "\n~~~~~~~~~~~~~~~~~~~~~~~~\n";

        // Write all initial probabilities into file by iterating over init_probs map:
        out << "\nInitial Probabilities:\n";
        
        for (auto e = init_probs.begin(); e != init_probs.end(); ++e) {
            out << get_state(e->first) << "\t" << e->second << "\n";  
        }
        out.close();
      }
      else {std::cerr << "Unable to open '" << filename << "'\n";}
     }
    else {std::cerr << "Memory error: Pointer to null. Could not read in file.\n";}
  }


  /**
    @brief Finds all bigrams that occurred in the file and counts them.
  */
  void get_bigrams () {   
    // The very first tag is in position[0].
    std::string previous_tag = pos_seq[0];      // Saves the tag that occurred before the current word 
    
    // Iterate over pos sequence, starting at the second word:
    for (unsigned i = 1; i < token_counter; ++i) {     
        
      // Pointer to the bigram that consists of the previous and current tag:
      auto bigram = std::make_tuple(previous_tag, pos_seq[i]);
      // Increment counter and make the current tag the new previous tag:
      ++bigram_counter_map[bigram];
      previous_tag = pos_seq[i];
    }
  }


  /**
    @brief Computes the maximum likelihood estimation for bigrams. If smoothing is set to true (default),
    it will compute probabilities for all possible bigrams using Add-One-Smoothing. If a number is specified, 
    it will use this number as smoothing summand (instead of one). If smoothing is set to false,
    it will only compute probabilities for the bigrams that occurred in the file.
  */
  void compute_mle() {
    // Smoothing = true, compute all possible bigram probabilities:
    if (smoothing == true) {
        
      // The vocabulary size multiplied with specified smoothing summand:
      float V = (float(type_counter_map.size()) * smoothing_summand);

      // Then iterate twice over pos_seq to get all possible bigram combinations:
      for (auto p = pos_seq.begin(); p < pos_seq.end(); ++p) {
        for (auto q = pos_seq.begin(); q < pos_seq.end(); ++q) {

          std::string tag1 = *p;                          // First tag of bigram
          std::string tag2 = *q;                          // Second tag of bigram
          BigramTuple bigram = std::make_pair(tag1,tag2); // Current bigram

          // Pointer to current bigram in bigram_counter_map:
          auto it = bigram_counter_map.find(bigram);
          
          // If we found an entry for the bigram, use its count for the computation:
          if (it != bigram_counter_map.end()) {
            ml_map[bigram] = (double(it->second) + smoothing_summand) / double(type_counter_map[tag1] + V);
          }
          // If we did not find the bigramm in the map, it did not occur in the training data:
          else {
            // Computation: The maximum likelihood estimation of a bigram consisting of
            // tag1 and tag2 is the the smoothing summand (=1 in Laplace Smoothing) divided by
            // the number of all occurrences of word 1 plus the size of the vocabulary.
            ml_map[bigram] = smoothing_summand / double(type_counter_map[tag1] + V);
          }
        }
      }
    }

    // Smoothing = false, compute the probability of all bigrams that occurred in the training data: 
    else {
      for (auto e = bigram_counter_map.begin(); e != bigram_counter_map.end(); ++e) {
        std::string tag1;       // First tag of bigram
        std::string tag2;       // Second tag of bigram
        
        // Use std::tie to unpack tuple into variable:
        std::tie (tag1, tag2) = (e->first);
        
        // Computation: The maximum likelihood estimation of a bigram consisting of
        // tag1 and tag2 is the number of all occurrences of this bigramm divided by
        // the number of all occurrences of word 1.
        ml_map[(e->first)] = double(e->second) / double(type_counter_map[tag1]);
      }
    }
  }


  /**
    @brief Computes the word likelihood for all observed words in the training data.
  */
  void compute_word_likelihood() {
      
    // Iterate over all observation tuples:
    for (auto e = observation_counter_map.begin(); e != observation_counter_map.end(); ++e) {
      std::string pos_tag;        // First tuple component: pos tag
      std::string word;           // Second tuple component: word
      
      // Use std::tie to unpack tuple into variable:
      std::tie (pos_tag, word) = (e->first);
      
      // Computation: The word likelihood of a pos tag and word tuple is the number of all occurrences
      // of this tuple divided by the number of all occurrences of the pos tag.
      word_likelihood_map[(e->first)] = double(e->second) / double(type_counter_map[pos_tag]);
    }
  }


  /**
    @brief Computes initial probablities for all existing states. Gets the initial probability for a state
           by looking for a transition from BOS to the state. 
    @note If smoothing is disabled, it only saves the real transitions from BOS to a state
     and the other states will have a zero as probability in the map (from initializing them). If smoothing 
          is enabled, we take all the smoothed BOS-state-transitions and save them.
  */      
  void get_init_probs() {
      
    // We iterate over the map that saves the transitions:
    for (auto e = ml_map.begin(); e != ml_map.end(); ++e) {
      std::string tag1;                       // First tag of bigram
      std::string tag2;                       // Second tag of bigram
      std::tie (tag1, tag2) = (e->first);
      
      // If the first tag is BOS, the second tag is a start state:
      if (tag1 == "<BOS>") {
        // Check whether second tag is in init probs map:
        auto it = init_probs.find(tag2);      // Pointer to the position of the specified tag in the map
        
        // If it is in there, overwrite the zero probability with the probability of the transition:
        if (it != init_probs.end()) {
          init_probs[tag2] = e->second;
        }
      }
    } 
  }


  /**
    @brief Initializes initial probability map by writing zero in it for all states.
  */
  void initialize_init_probs() {
    for (auto e = pos_seq.begin(); e != pos_seq.end(); ++e) {
      init_probs.insert(std::make_pair((*e),0));
    }
  }


  /**
    @brief Creates the map that maps all pos tags to their state index.
  */
  void build_tag_state_map() {
    unsigned state_index = 0;     // Current state index  
    
    // Iterate over all pos-tags: 
    for (auto e = type_counter_map.begin(); e != type_counter_map.end(); ++e) {
      // Save a mapping from current tag to current index, then increment index:
      tag_state_map[e->first] = state_index;
      ++state_index;
    }
  }


  /**
    @brief Returns the state representation for a given pos-tag.
    @param tag The pos-tag we want to get the state for 
    @return The state of the tag
  */
  unsigned get_state(const std::string& tag) const 
  {
    // Pointer to the position of the specified tag in the map:
    auto it = tag_state_map.find(tag);
    
    if (it != tag_state_map.end()) return it->second;  
    else {
      std::cerr << "ERROR: No state representation for pos tag " << tag << " found.\n";
      exit(1);
    }
  }


  /**
    @brief Removes all duplicates from a string vector.
    @param vec The vector to remove duplicates from 
  */
  void remove_duplicates(StringVector& vec) {
    std::sort(vec.begin(), vec.end());                         // Sort elements 
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); // Remove duplicate elements
  }

};

#endif