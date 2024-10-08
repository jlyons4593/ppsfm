#include "dataStructures.hpp"
#include "Helper.hpp"
#pragma once
class EstimatedRobustViews{
  private:


DataStructures::SfMData data;
DataStructures::ComputedCameraPoints camera_variables;

Helper::InlierResults find_inliers(Eigen::MatrixXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points);

std::pair<Eigen::VectorXd, Eigen::RowVectorXd> compute_reproj(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points);

double computeScore(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points, std::vector<int> inliers);

  public:
    EstimatedRobustViews(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables,Eigen::VectorXi known_points,int new_view, int num_rejected, int level);

    Eigen::VectorXd best_estimate;
    Eigen::VectorXi best_inliers; 

};
