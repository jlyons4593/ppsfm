
#include "dataStructures.hpp"

#include "Helper.hpp"

#pragma once
class EstimatedRobustPoints{
  private:

  public:


DataStructures::SfMData data;
DataStructures::ComputedCameraPoints camera_variables;

Helper::InlierResults find_inliers(Eigen::MatrixXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points);

std::pair<Eigen::VectorXd, Eigen::RowVectorXd> compute_reproj(Eigen::VectorXd& estimation,Eigen::VectorXi& idx_view, Eigen::VectorXi& known_points);

double computeScore(Eigen::VectorXd& estimation,Eigen::MatrixXd& cams , Eigen::VectorXi& idx_view, Eigen::VectorXi& known_views, std::vector<int> inliers);

  public:
    EstimatedRobustPoints(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables,Eigen::VectorXi& known_points,int new_point, int num_rejected, int level);

    Eigen::VectorXd best_estimate;
    Eigen::VectorXi best_inliers; 
};
