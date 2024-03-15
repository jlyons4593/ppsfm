#include "dataStructures.hpp"
#pragma once
struct InlierResults {
    std::vector<int> inliers;
    double score;
};
class EstimatedRobustViews{

  private:

DataStructures::SfMData data;
DataStructures::ComputedCameraPoints camera_variables;

InlierResults find_inliers(Eigen::MatrixXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points);
// InlierResults find_inliers(const Eigen::MatrixXd& estimation, 
                           // const Eigen::MatrixXd& points, 
                           // const Eigen::VectorXi& idx_points, 
                           // const Eigen::VectorXi& idx_view,
                           // const Eigen::MatrixXd& normalisations, 
                           // const Eigen::MatrixXd& img_meas, 
                           // double threshold, 
                           // bool sort_inliers = false) ;

// Eigen::VectorXd compute_reproj(const Eigen::MatrixXd& estimation, 
//                                const Eigen::MatrixXd& points, 
//                                const Eigen::VectorXi& idx_points, 
//                                const Eigen::VectorXi& idx_view, 
//                                const Eigen::MatrixXd& normalisations, 
//                                const Eigen::MatrixXd& img_meas);

std::pair<Eigen::VectorXd, Eigen::RowVectorXd> compute_reproj(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points);
  public:
    EstimatedRobustViews(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables,Eigen::VectorXi known_points,int new_view, int num_rejected, int level);


};
