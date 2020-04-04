/*
 * Alevin Efficient Data Storage (EDS) reader
 *
 * Author: Avi Srivastava
 * Last modified: August 13, 2019
 * License: LGPL (>= 3)
 *
 */

#include <Rcpp.h>
#include <zlib.h>

using namespace Rcpp;

// C++ internal function to figure out the spaces to reserve
size_t getReserveSpaces(size_t numOfGenes, size_t numOfOriginalCells,
                      Rcpp::IntegerVector& bitVecLengths,
                      std::string& countMatFilename,
                      bool tierImport = false) {

  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(countMatFilename.c_str(), "rb") ;
  
  // We are storing the bit vector in u8 so total number of u8 = numGenes/8
  size_t numFlags = std::ceil(numOfGenes / 8.0);

  // vector for storing the bitvector flags
  std::vector<uint8_t> alphasFlag (numFlags, 0);

  // getting the sizs of u8 and float 32
  size_t totalSpace { 0 };
  size_t flagSize = sizeof(decltype(alphasFlag)::value_type);
  size_t elSize = sizeof(float);
  if(tierImport) { elSize = sizeof(uint8_t); }

  // iterating over cells
  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {
    // reading bitvectors
    gzread(fileHandler, reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);
    size_t numOfExpGenes { 0 };

    for (size_t j = 0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i = 0; i < 8; i++){
        // counting positions only if the flag is set
        if (flag & (128 >> i)) {
          numOfExpGenes += 1;
        }
      }
    }

    // skipping the expression values and saving the counts for numOfExpGenes
    gzseek(fileHandler, elSize * numOfExpGenes, SEEK_CUR);
    bitVecLengths[ cellId + 1 ] = ( numOfExpGenes + bitVecLengths[ cellId ] );
    totalSpace += numOfExpGenes;
  }

  return totalSpace;
}

template <typename T>
bool populateCounts(size_t elSize, size_t numExpGenes, gzFile& fileHandler,
                    size_t& valCounter, size_t totalSpace, 
                    Rcpp::NumericVector& values) {
    // reading in the expression
    std::vector<T> alphasSparse(numExpGenes);
    gzread(fileHandler, reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

    // saving the positions and expression
    for (size_t i = 0; i < numExpGenes; i++) {
      if ( valCounter >= totalSpace ) {
        return false;
      }

      values[valCounter] = alphasSparse[i];
      valCounter += 1;
    }

    return true;
}

// [[Rcpp::export]]
SEXP getSparseMatrix(size_t numOfGenes, size_t numOfOriginalCells, 
                     std::string countMatFilename,
                     bool tierImport = false) {
  Rcpp::S4 mat("dgCMatrix");
  
  // initializing vector to store bitvecSpaces
  Rcpp::IntegerVector bitVecLengths(numOfOriginalCells + 1, 0);
  size_t totalSpace = getReserveSpaces(numOfGenes, numOfOriginalCells, 
                                       bitVecLengths, countMatFilename,
                                       tierImport);

  // initializing sparse matrix
  typedef Rcpp::NumericVector ValuesT;
  ValuesT values(totalSpace, 0.0);

  // initializing sparse matrix indices
  typedef Rcpp::IntegerVector IndicesT;
  IndicesT indices(totalSpace, 0);

  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(countMatFilename.c_str(), "rb") ;
  
  // We are storing the bit vector in u8 so total number of u8 = numGenes/8
  size_t numFlags = std::ceil(numOfGenes / 8.0);

  // vector for storing the bitvector flags
  std::vector<uint8_t> alphasFlag (numFlags, 0);

  // getting the sizs of u8 and float 32
  size_t flagSize = sizeof(decltype(alphasFlag)::value_type);
  size_t elSize = sizeof(float);
  if(tierImport) { elSize = sizeof(uint8_t); }

  size_t valCounter { 0 };
  // iterating over cells
  for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId) {
    // reading bitvectors
    gzread(fileHandler, reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);

    // iterating over u8 flags for bitvectors
    size_t numExpGenes { 0 };
    for (size_t j = 0; j < alphasFlag.size(); j++) {
      uint8_t flag = alphasFlag[j];

      for (size_t i = 0; i < 8; i++){
        // extracting positions only if the flag is set
        if (flag & (128 >> i)) {
          if ( valCounter + numExpGenes >= totalSpace ) { 
            return Rcpp::List(); 
          }
          
          size_t offset = i + (8 * j);
          indices[ valCounter + numExpGenes ] = offset;
          numExpGenes += 1;
        }
      }
    }
    
    bool parseOk = tierImport ? populateCounts<uint8_t>(elSize, numExpGenes, fileHandler,
                                                        valCounter, totalSpace, values) :
                                populateCounts<float>(elSize, numExpGenes, fileHandler,
                                                      valCounter, totalSpace, values);

    if ( not parseOk ) {
        return Rcpp::List();
    }
  }

  // code in-parts taken from https://github.com/LTLA/beachmat/blob/master/inst/include/beachmat/output/Csparse_writer.h#L268
  mat.slot("Dim") = Rcpp::IntegerVector::create(numOfGenes, numOfOriginalCells);

  // Setting p
  mat.slot("p") = bitVecLengths;

  // Setting 'x'.
  mat.slot("x") = values;

  // Setting 'i'.
  mat.slot("i") = indices;

  return SEXP(mat);
}
