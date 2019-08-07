/*
 * Alevin Efficient Data Storage (EDS) reader
 * 
 * Author: Avi Srivastava
 * Last modified: August 7, 2019
 * License: LGPL (>= 3)
 *
 */

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <zlib.h>
using namespace Rcpp;

// C++ internal function to figure out the spaces to reserve
void getReserveSpaces(size_t numOfGenes, size_t numOfOriginalCells,
                      std::vector<size_t>& bitVecLngths,
                      std::string& countMatFilename) {

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
    size_t numOfExpGenes { 0 };

    for (size_t j=0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i=0; i<8; i++){
        // counting positions only if the flag is set
        if (flag & (128 >> i)) {
          numOfExpGenes += 1;
        }
      }
    }

    // skipping the expression values and saving the counts for numOfExpGenes
    gzseek(fileHandler, elSize * numOfExpGenes, SEEK_CUR);
    bitVecLngths[cellId] = numOfExpGenes;
  }
}

// [[Rcpp::export]]
SEXP getSparseMatrix(size_t numOfGenes, size_t numOfOriginalCells, std::string countMatFilename) {
  // initializing vector to store bitvecSpaces
  std::vector<size_t> reserveSpaces(numOfOriginalCells);
  getReserveSpaces( numOfGenes, numOfOriginalCells, reserveSpaces, countMatFilename );

  // initializing sparse matrix
  typedef Eigen::SparseMatrix<float> SpMatT;
  SpMatT spMat(numOfGenes, numOfOriginalCells);
  spMat.reserve(reserveSpaces);

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
    std::vector<uint32_t> indices;
    for (size_t j=0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i=0; i<8; i++){
        // extracting positions only if the flag is set
        if (flag & (128 >> i)) {
          indices.emplace_back( i+(8*j) );
        }
      }
    }
    
    // reading in the expression
    size_t numExpGenes { indices.size() };
    std::vector<float> alphasSparse(numExpGenes);
    gzread(fileHandler, reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

    // saving the positions and expression
    for (size_t i = 0; i<numExpGenes; i++) {
      spMat.insert( indices[i], cellId ) = alphasSparse[i];
    }
  }

  S4 Xout(wrap(spMat));
  return Xout;
}
