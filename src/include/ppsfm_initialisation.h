#include "dataStructures.hpp"
#include "estimatedFundamentalMatrix.h"
#include "logger.h"
#include <array>
#include "options.h"
#include <Eigen/Dense>
#include <iostream>
#include <ostream>
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

    DataStructures::ComputedCameraPoints compute_cams(std::vector<int> initial_views, std::vector<int> initial_points, Eigen::MatrixXd fundamental_matrix){
        DataStructures::ComputedCameraPoints camera_points; 
        Eigen::Vector3i first_view(3*initial_views[0],3*initial_views[0]+1, 3*initial_views[0]+2); 
        Eigen::Vector3i second_view(3*initial_views[1],3*initial_views[1]+1, 3*initial_views[1]+2); 

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(fundamental_matrix,Eigen::ComputeFullV);
        Eigen::MatrixXd V = svd.matrixV();

        // Epipole is the 3rd column
        Eigen::Vector3d epipole = V.col(2);
        
        Eigen::MatrixXd proj_depths = Eigen::MatrixXd::Ones(2, initial_points.size());

        for(int point_idx = 0; point_idx<initial_points.size(); point_idx++){
            Eigen::Vector3d vector_from_measurement;
            for(int i = 0; i < 3; ++i) {
                vector_from_measurement(i) = measurement(second_view[i], initial_points[point_idx]);
            }

            // Compute the cross product
            Eigen::VectorXd point_epipole = vector_from_measurement.cross(epipole);
            int col_idx = initial_points[point_idx]; 
            Eigen::Vector3d meas_vec;
            for (int i = 0; i < 3; ++i) {
                meas_vec(i) = measurement(first_view(i), col_idx);
            }

            double numerator = (meas_vec.transpose() * fundamental_matrix * point_epipole)(0,0); // This extracts the scalar value
            double denominator = (point_epipole.transpose() * point_epipole)(0, 0); // Also a scalar

            double value = std::abs(numerator / denominator);
            proj_depths(1,point_idx)= value;
        }
        double sumOfFirstTwo = proj_depths.row(1).segment(0, 2).sum();

        // Normalize the second row based on this sum and scale by 2
        proj_depths.row(1) = (proj_depths.row(1) / sumOfFirstTwo) * 2;
        for (int i = 2; i < proj_depths.cols(); ++i) {
            // Compute the sum of the first two elements in the column
            double sum = proj_depths.block(0, i, 2, 1).sum();

            // Normalize the column by the sum and multiply by 2
            if (sum != 0) { // preventing division by zero
                proj_depths.col(i) = proj_depths.col(i) / sum * 2;
            }
        }
        pathway.push_back(initial_points[0]);
        pathway.push_back(-initial_views[0]);

        pathway.push_back(initial_points[1]);
        pathway.push_back(-initial_views[1]);

        for(size_t i = 2; i < initial_points.size(); ++i) {
            pathway.push_back(initial_points[i]);
        }
        camera_points.fixed.resize(initial_points.size() + 2);
        camera_points.fixed[0] = initial_views[0];
        camera_points.fixed[1] = initial_points[0];
        camera_points.fixed[2] = initial_views[1]; 
        std::vector<int> combinedPoints = {initial_points[0], initial_points[1]}; 
        camera_points.fixed[3] = combinedPoints;

        for (Eigen::Index i = 4; i < camera_points.fixed.size(); ++i) {
            camera_points.fixed[i] = initial_views; 
        }
        std::vector<int> combined_indices;
        for(auto i: first_view){
            combined_indices.push_back(i);
        }
        combined_indices.insert(combined_indices.end(), second_view.begin(), second_view.end());

        // Create a scaled_measurement matrix with the correct dimensions
        Eigen::MatrixXd scaled_measurement(combined_indices.size(), initial_points.size());

        // Fill in scaled_measurement by manually selecting the appropriate elements
        for (size_t i = 0; i < combined_indices.size(); ++i) {
            for (size_t j = 0; j < initial_points.size(); ++j) {
                scaled_measurement(i, j) = measurement(combined_indices[i], initial_points[j]);
            }
        }
        assert(proj_depths.rows() == 2);
        //
        Eigen::MatrixXd adjusted_proj_depths(6, proj_depths.cols());
        for (int i = 0; i < proj_depths.cols(); ++i) {
            for (int j = 0; j < 3; ++j) { // First 3 rows are copies of the first row of proj_depths
                adjusted_proj_depths.row(j)(i) = proj_depths.row(0)(i);
            }
            for (int j = 3; j < 6; ++j) { // Next 3 rows are copies of the second row of proj_depths
                adjusted_proj_depths.row(j)(i) = proj_depths.row(1)(i);
            }
        }
        //
        // Perform element-wise multiplication
        scaled_measurement = scaled_measurement.array() * adjusted_proj_depths.array();
        Eigen::JacobiSVD<Eigen::MatrixXd> scaled_svd(scaled_measurement, Eigen::ComputeFullV | Eigen::ComputeThinU);

        // Extract the matrices U, S, and V
        Eigen::MatrixXd scaled_U = scaled_svd.matrixU();
        Eigen::MatrixXd scaled_S = scaled_svd.singularValues().asDiagonal();
        Eigen::MatrixXd scaled_V = scaled_svd.matrixV();


        Eigen::VectorXd sqrt_singularValues = scaled_svd.singularValues().array().sqrt();

        Eigen::MatrixXd sqrt_S_4 = sqrt_singularValues.head(4).asDiagonal();
        

        // Now perform the multiplication with the first four columns of U
        camera_points.cameras = scaled_U.leftCols(4) * sqrt_S_4;
        camera_points.points = (sqrt_S_4 * scaled_V.leftCols(4).transpose());       

        // HACKY CONVERSION PLEASE FIX
        Eigen::Map<Eigen::VectorXi> temp_vector(pathway.data(), pathway.size());
        camera_points.pathway = temp_vector;
        return camera_points;

    }
public:

    Initialisation(const Eigen::MatrixXd& measurement, const Eigen::MatrixXd& visible, 
                   const Eigen::MatrixXd& view_pairs, const Eigen::VectorXd& affinity)
        : measurement(measurement), visible(visible), view_pairs(view_pairs),
        affinity(affinity) {}
    
    DataStructures::ComputedCameraPoints getCameraPoints(){
        return camera_variables;
    }

    void process();
    // void process() {
    
};
