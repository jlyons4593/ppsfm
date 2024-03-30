
#include <Eigen/Dense>

#pragma once
namespace Helper{
struct InlierResults {
    std::vector<int> inliers;
    double score;
    Eigen::VectorXd reproj_errors;
};


  Eigen::VectorXi random_subset(const Eigen::VectorXi& complete_set, int num_rejected, int num_sample);

} 

