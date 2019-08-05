/*
 * Alevin Efficient Data Storage (EDS) reader
 * 
 * Author: Avi Srivastava
 * Last modified: August 1, 2019
 * License: LGPL (>= 3)
 *
 */

#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::vector<size_t>> getPositions(size_t numOfGenes,
					      size_t numOfOriginalCells, 
					      std::string countMatFilename) {

  std::vector<std::vector<size_t>> positions;
  std::ifstream fileHandler(countMatFilename.c_str(), std::ios::in | std::ios::binary);
  
  size_t numFlags = std::ceil(numOfGenes / 8.0);
  std::vector<uint8_t> alphasFlag (numFlags, 0);
  size_t flagSize = sizeof(decltype(alphasFlag)::value_type);
  size_t elSize = sizeof(float);

  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {
    fileHandler.read(reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);

    std::vector<size_t> indices;
    for (size_t j=0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i=0; i<8; i++){
        if (flag & (128 >> i)) {
          indices.emplace_back( i+(8*j) );
        }
      }
    }

    fileHandler.ignore(elSize * indices.size());
    positions.emplace_back(indices);
  }

  return positions;
}

// [[Rcpp::export]]
std::vector<std::vector<float>> getExpression(size_t numOfGenes,
					      size_t numOfOriginalCells, 
                                              std::vector<std::vector<size_t>>& positions,
                                              std::string countMatFilename) {
  std::vector<std::vector<float>> counts;
  std::ifstream fileHandler(countMatFilename.c_str(), std::ios::in | std::ios::binary) ;

  size_t numFlags = std::ceil(numOfGenes / 8.0);
  std::vector<float> alphasSparse;
  alphasSparse.reserve(numFlags/2);
  size_t flagSize = sizeof(uint8_t);
  size_t elSize = sizeof(decltype(alphasSparse)::value_type);

  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {

    size_t numExpGenes { positions[cellId].size() };
    fileHandler.ignore(flagSize * numFlags);
    alphasSparse.clear();
    alphasSparse.resize(numExpGenes);
    fileHandler.read(reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

    counts.emplace_back(alphasSparse);
  }

  return counts;
}
