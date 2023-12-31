#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mat.h"
// void matread(const char *file, std::vector<double>& v)
// {
//     // open MAT-file
//     MATFile *pmat = matOpen(file, "r");
//     if (pmat == NULL) return;
//
//     // extract the specified variable
//     mxArray *arr = matGetVariable(pmat, "LocalDouble");
//     if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr)) {
//         // copy data
//         mwSize num = mxGetNumberOfElements(arr);
//         double *pr = mxGetPr(arr);
//         if (pr != NULL) {
//             v.reserve(num); //is faster than resize :-)
//             v.assign(pr, pr+num);
//         }
//     }
//
//     // cleanup
//     mxDestroyArray(arr);
//     matClose(pmat);
// }
void matread(const char *file, Eigen::VectorXd& eigenVec)
{
    // Open MAT-file
    MATFile *pmat = matOpen(file, "r");
    if (pmat == NULL) return;

    // Extract the specified variable
    mxArray *arr = matGetVariable(pmat, "LocalDouble");
    if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr)) {
        // Copy data
        mwSize num = mxGetNumberOfElements(arr);
        double *pr = mxGetPr(arr);
        if (pr != NULL) {
            // Use Eigen::Map to map the data directly into an Eigen vector
            eigenVec = Eigen::Map<Eigen::VectorXd>(pr, num);
        }
    }

    // Cleanup
    mxDestroyArray(arr);
    matClose(pmat);
}

void pipeline(Eigen::SparseMatrix<double> measurements, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>){
    
}

int main()
{ 
    Eigen::VectorXd eigenVector;
    matread("Data/data.mat", eigenVector);
    for (size_t i=0; i<eigenVector.size(); ++i)
        std::cout << eigenVector[i] << std::endl;
    
    std::cout<< "Hello World"<< std::endl;
    return 0;
}
