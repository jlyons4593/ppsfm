#include <iostream>
#include "mat.h"
void matread(const char *file, std::vector<double>& v)
{
    // open MAT-file
    MATFile *pmat = matOpen(file, "r");
    if (pmat == NULL) return;

    // extract the specified variable
    mxArray *arr = matGetVariable(pmat, "LocalDouble");
    if (arr != NULL && mxIsDouble(arr) && !mxIsEmpty(arr)) {
        // copy data
        mwSize num = mxGetNumberOfElements(arr);
        double *pr = mxGetPr(arr);
        if (pr != NULL) {
            v.reserve(num); //is faster than resize :-)
            v.assign(pr, pr+num);
        }
    }

    // cleanup
    mxDestroyArray(arr);
    matClose(pmat);
}

int main()
{
    std::vector<double> v;
    matread("Data/data.mat", v);
    for (size_t i=0; i<v.size(); ++i)
        std::cout << v[i] << std::endl;
    std::cout<< "hello world!"<< std::endl;
    return 0;
}
// Sample Eigen Code for using simple matrices
// #include <iostream>
// #include <Eigen/Dense>

// int main() {

    // Eigen::Matrix3f matrix_3x3;  // 3x3 matrix with float elements matrix_3x3 << 1.0f, 2.0f, 3.0f,
                  // 4.0f, 5.0f, 6.0f,
                  // 7.0f, 8.0f, 9.0f;
// 
    // Eigen::Vector3f vec;  // 3D vector with float elements
    // vec << 1.0f, 2.0f, 3.0f;
// 
    // Eigen::Vector3f result = matrix_3x3 * vec;  // Result is also a float vector
// 
    // std::cout << "Matrix:\n" << matrix_3x3 << "\n";
    // std::cout << "Vector:\n" << vec << "\n";
    // std::cout << "Result:\n" << result << "\n";
// 
    // return 0;
// }
