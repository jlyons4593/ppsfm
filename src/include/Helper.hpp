
#include <Eigen/Dense>

#pragma once
namespace Helper{


  Eigen::VectorXi random_subset(const Eigen::VectorXi& complete_set, int num_rejected, int num_sample);

} 

