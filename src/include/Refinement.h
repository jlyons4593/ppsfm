#include "dataStructures.hpp"
#pragma once
class Refinement{
  private:
    DataStructures::SfMData data;
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> visible;
    DataStructures::ComputedCameraPoints camera_variables;
    void reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points);
    void reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras);
  public:
    Eigen::MatrixXd getCameras();
    Eigen::MatrixXd getPoints();
    
    Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> inliers, int init_refine, int last_path, bool start_cameras, int type);


};
