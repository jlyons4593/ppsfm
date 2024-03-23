#include "dataStructures.hpp"
#pragma once
class Refinement{
  private:
    Eigen::Array<std::variant<int, std::vector<int>, Eigen::VectorXi>, Eigen::Dynamic, 1> fixed; 
    DataStructures::SfMData data;
    DataStructures::ComputedCameraPoints camera_variables;
    void reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras);
    Eigen::MatrixXd reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points);

  public:
    Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables , int init_refine, int last_path, bool start_cameras, int type);


};
