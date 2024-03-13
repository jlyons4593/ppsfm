#include "dataStructures.hpp"
#include "estimatedFundamentalMatrix.h"
#include <Eigen/Dense>
#include <vector>

// Assuming Options is a struct that contains all necessary fields

#pragma once
class Initialisation{
private:

    DataStructures::ComputedCameraPoints camera_variables;
    Eigen::MatrixXd measurement; // 3FxN matrix
    Eigen::MatrixXd visible;     // FxN binary mask matrix
    Eigen::MatrixXd view_pairs;  // Kx2 matrix
    Eigen::VectorXd affinity;    // Kx1 vector
    // Assuming estimated_views is a vector of integers representing the IDs of views that have been estimated
    Eigen::VectorXi estimated_views;

    Eigen::MatrixXd cameras; // 6x4 matrix for projective estimation of initial cameras
    Eigen::MatrixXd points;  // 4xK matrix for projective estimation of initial points
    std::vector<int> pathway; // Order in which views and points are added
    // Assuming 'fixed' is a vector of vectors, representing the points or views used in constraints

    DataStructures::ComputedCameraPoints compute_cams(std::vector<int> initial_views, std::vector<int> initial_points, Eigen::MatrixXd fundamental_matrix);
    
public:

    Initialisation(const Eigen::MatrixXd& measurement, const Eigen::MatrixXd& visible, 
                   const Eigen::MatrixXd& view_pairs, const Eigen::VectorXd& affinity);
    
    DataStructures::ComputedCameraPoints getCameraPoints();

    void process();
    
};
