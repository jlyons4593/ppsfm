#include "dataStructures.hpp"

#pragma once
class EstimatedPoints{
  private:
  public:

    Eigen::MatrixXd sys;  

    // Initialize linear equality constraint (fixes sum of the fixed projective depths)
    Eigen::RowVectorXd con; 

    // Initialize positivity constraint (projective depths)
    Eigen::MatrixXd pos; 

    Eigen::VectorXd estim;

    EstimatedPoints(DataStructures::SfMData& data,Eigen::MatrixXd& cameras, Eigen::VectorXi& idx_views,int new_point, int num_fixed);
};
