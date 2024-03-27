#include "diagnosis.h"
#include "logger.h"
Diagnosis::Diagnosis(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& inliers, DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables, int iteration){
    int fig_base = 0;
    int num_views = data.visible.rows();
    int num_points = data.visible.cols();
    std::cout<<num_views<< " "<<num_points<<std::endl;

    auto& valid_cams = camera_variables.negative_pathway;
    std::cout<<valid_cams<<std::endl;
    auto& valid_points = camera_variables.positive_pathway;
    
    Logger::logSubsection("Diagnosis");

    Eigen::MatrixXd unknowns = data.visible;

    // Set unknowns to 0 for valid cameras and valid points

    // Set unknowns to false for valid cameras and valid points
    for (int i = 0; i < valid_cams.size(); ++i) {
        for (int j = 0; j < valid_points.size(); ++j) {
            unknowns(valid_cams(i), valid_points(j)) = 0;
        }
    }
    // Update inliers by removing unknowns
    inliers = (inliers.array() != 0 && unknowns.array() == 0);

    // Count the number of inliers
    int num_inliers = (inliers.array() != 0).cast<int>().sum();

    // Count the number of visible points that are not unknowns
    Eigen::MatrixXd visible_not_unknowns = (data.visible.array() != 0 && unknowns.array() == 0).cast<double>();
    int num_visible_not_unknowns = (visible_not_unknowns.array() != 0).cast<int>().sum();

    // Calculate the inlier ratio
    double inlier_ratio = static_cast<double>(num_inliers) / num_visible_not_unknowns * 100;

    double view_percentage =(double(valid_cams.size())/double(num_views))*100;
    double point_percentage =(double(valid_points.size())/double(num_points))*100;
    std::cout<<"Estimated "<<valid_cams.size()<<"/"<<num_views<<" views ("<<view_percentage<<"%) and "<<valid_points.size()<<"/"<<num_points<<" ("<<point_percentage<< "%) with "<<inlier_ratio<<"% inliers."<<std::endl;

}
