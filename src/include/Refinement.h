#include "dataStructures.hpp"
#pragma once
class Refinement{
  private:
    DataStructures::SfMData data;
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> visible;
    DataStructures::ComputedCameraPoints camera_variables;
void reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points, std::vector<std::vector<int>> fixed);

void reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras, std::vector<std::vector<int>> fixed);
    Eigen::VectorXi pathway_segment;
    std::vector<std::vector<int>> fixed;
  public:
    Eigen::MatrixXd getCameras();
    Eigen::MatrixXd getPoints();
    
    Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> visible,Eigen::VectorXi pathway_segment,std::vector<std::vector<int>> new_fixed,  bool start_cameras, int type);


};
