#include "dataStructures.hpp"
#include "mat.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

#pragma once
namespace DataReader {

// Overload the << operator to print the Matrices struct
inline std::ostream &operator<<(std::ostream &os,
                                const DataStructures::InputMatrices &matrices) {
  os << "Measurements:\n" << matrices.measurements << "\n\n";
  os << "Image Size:\n" << matrices.image_size << "\n\n";
  os << "Centers:\n" << matrices.centers << "\n";
  return os;
}

inline DataStructures::InputMatrices matread(const char *file) {
  DataStructures::InputMatrices matrices;

  // Open MAT-file
  MATFile *pmat = matOpen(file, "r");
  if (pmat == NULL) {
    std::cerr << "issue loading file" << std::endl;
    return matrices;
  }
  // Helper lambda to load a sparse matrix
  auto loadSparseMatrix = [&](const char *name,
                              Eigen::SparseMatrix<double> &matrix) {
    mxArray *arr = matGetVariable(pmat, name);
    if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr) && mxIsSparse(arr)) {
      // Get the sparse matrix components
      mwSize numRows = mxGetM(arr);
      mwSize numCols = mxGetN(arr);
      double *pr = mxGetPr(arr);
      mwIndex *ir = mxGetIr(arr); // Row indices
      mwIndex *jc = mxGetJc(arr); // Column index pointers

      // Reserve space for non-zeros
      matrix.resize(numRows, numCols);
      matrix.reserve(mxGetNzmax(arr)); // Reserve space for non-zero elements

      // Fill the sparse matrix
      for (mwSize col = 0; col < numCols; ++col) {
        for (mwIndex index = jc[col]; index < jc[col + 1]; ++index) {
          matrix.insert(ir[index], col) = pr[index];
        }
      }
    }
    mxDestroyArray(arr);
  };

  // Load the sparse matrix measurements
  loadSparseMatrix("measurements", matrices.measurements);

  // Helper lambda to load a matrix
  auto loadMatrix = [&](const char *name, Eigen::MatrixXd &matrix) {
    mxArray *arr = matGetVariable(pmat, name);
    if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr)) {
      mwSize numElements = mxGetNumberOfElements(arr);
      mwSize numRows = mxGetM(arr);
      double *pr = mxGetPr(arr);
      if (pr != NULL) {
        // Map the data directly into an Eigen matrix
        matrix =
            Eigen::Map<Eigen::MatrixXd>(pr, numRows, numElements / numRows);
      }
      mxDestroyArray(arr);
    }
  };

  // Load each matrix
  loadMatrix("image_size", matrices.image_size);
  loadMatrix("centers", matrices.centers);

  // Cleanup
  matClose(pmat);

  return matrices;
}
} // namespace DataReader
