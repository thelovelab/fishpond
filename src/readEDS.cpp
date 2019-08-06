/*
 * Alevin Efficient Data Storage (EDS) reader
 * 
 * Author: Avi Srivastava
 * Last modified: August 6, 2019
 * License: LGPL (>= 3)
 *
 */

#include <Rcpp.h>
#include <zlib.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::vector<size_t>> getPositions(size_t numOfGenes, size_t numOfOriginalCells, 
                                               std::string countMatFilename) {
  // initializing vector to store positions
  std::vector<std::vector<size_t>> positions;

  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(countMatFilename.c_str(), "rb") ;
  
  // We are storing the bit vector in u8 so total number of u8 = numGenes/8
  size_t numFlags = std::ceil(numOfGenes / 8.0);

  // vector for storing the bitvector flags
  std::vector<uint8_t> alphasFlag (numFlags, 0);

  // getting the sizs of u8 and float 32
  size_t flagSize = sizeof(decltype(alphasFlag)::value_type);
  size_t elSize = sizeof(float);

  // iterating over cells
  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {
    // reading bitvectors
    gzread(fileHandler, reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);

    // iterating over u8 flags for bitvectors
    std::vector<size_t> indices;
    for (size_t j=0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i=0; i<8; i++){
        // extracting positions only if the flag is set
        if (flag & (128 >> i)) {
          indices.emplace_back( i+(8*j) );
        }
      }
    }

    // skipping the expression values and saving the positions
    gzseek(fileHandler, elSize * indices.size(), SEEK_CUR);
    positions.emplace_back(indices);
  }

  return positions;
}

// [[Rcpp::export]]
std::vector<std::vector<float>> getExpression(size_t numOfGenes, size_t numOfOriginalCells, 
                                              std::vector<std::vector<size_t>>& positions,
                                              std::string countMatFilename) {
  // initializing vector to store expression counts
  std::vector<std::vector<float>> counts;

  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(countMatFilename.c_str(), "rb") ;

  // We are storing the bit vector in u8 so total number of u8 = numGenes/8
  size_t numFlags = std::ceil(numOfGenes / 8.0);
  
  // vector for storing the per cell expresssion
  std::vector<float> alphasSparse;
  alphasSparse.reserve(numFlags/2);

  // getting the sizs of u8 and float 32
  size_t flagSize = sizeof(uint8_t);
  size_t elSize = sizeof(decltype(alphasSparse)::value_type);

  // iterating over u8 flags for bitvectors
  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {

    // skipping the bit vectors
    size_t numExpGenes { positions[cellId].size() };
    gzseek(fileHandler, flagSize * numFlags, SEEK_CUR);
    
    // reading in the expression
    alphasSparse.clear();
    alphasSparse.resize(numExpGenes);
    gzread(fileHandler, reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

    // saving the expression
    counts.emplace_back(alphasSparse);
  }

  return counts;
}
