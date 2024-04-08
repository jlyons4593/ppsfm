
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#pragma once
namespace DataStructures {

struct ModelData{
    std::vector<double> timings;
    Eigen::MatrixXd initial_cameras;
    Eigen::MatrixXd initial_points;
    Eigen::VectorXi initial_pathway;
    std::vector<std::vector<int>> initial_fixed;
    Eigen::MatrixXd cameras;
    Eigen::MatrixXd points;
    Eigen::VectorXi pathway;
    std::vector<std::vector<int>> fixed;
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> inliers; 
};
struct InputMatrices {
  Eigen::SparseMatrix<double> measurements;
  Eigen::MatrixXd image_size;
  Eigen::MatrixXd centers;
};
struct ComputedCameraPoints {
  Eigen::MatrixXd cameras;
  Eigen::MatrixXd points;
  Eigen::RowVectorXi pathway;
  Eigen::RowVectorXi positive_pathway;
  Eigen::RowVectorXi negative_pathway;

  // Find a way to remove the std::vector from this
  std::vector<std::vector<int>> fixed;
};

inline void printColsRows(Eigen::MatrixXd matrix, std::string matrix_name) {

  std::cout << matrix_name + " rows: " << matrix.rows() << std::endl
            << matrix_name + " columns: " << matrix.cols() << std::endl;
}
struct ViewpairAffinity {
  Eigen::MatrixXd view_pairs;
  Eigen::VectorXd Affinity;
  ViewpairAffinity() {
    view_pairs = Eigen::MatrixXd(0, 0);
    Affinity = Eigen::VectorXd(0);      
  }
};


struct FinalData {
  Eigen::MatrixXd centers;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> visible; 
  Eigen::MatrixXd cost_function_data;    
  Eigen::MatrixXd normalised_measurements;
  Eigen::MatrixXd normalisations; 
  Eigen::MatrixXd
      image_measurements; 
  Eigen::MatrixXd
      pseudo_inverse_measurements; 
  Eigen::Array<bool, 1, Eigen::Dynamic>
      removed_points; 

  FinalData() = default;

};
struct SfMData {
  Eigen::MatrixXd centers;
  Eigen::MatrixXd visible; 
  Eigen::MatrixXd cost_function_data;    
  Eigen::MatrixXd normalised_measurements;
  Eigen::MatrixXd normalisations; 
  Eigen::MatrixXd
      image_measurements; 
  Eigen::MatrixXd
      pseudo_inverse_measurements; 
  Eigen::Array<bool, 1, Eigen::Dynamic>
      removed_points; 

  SfMData() = default;

  SfMData(int F, int N)
      : visible(Eigen::MatrixXd::Zero(F, N)),
        cost_function_data(Eigen::MatrixXd::Zero(
            F, N)), 
        normalised_measurements(Eigen::MatrixXd::Zero(3 * F, N)),
        normalisations(Eigen::MatrixXd::Zero(3 * F, 3)),
        image_measurements(Eigen::MatrixXd::Zero(3 * F, N)),
        pseudo_inverse_measurements(Eigen::MatrixXd::Zero(F, 3 * N)),
        removed_points(Eigen::Array<bool, 1, Eigen::Dynamic>()) {}
};

struct Model {
  Eigen::MatrixXi inliers;             // FxN inliers matrix
  Eigen::MatrixXd cameras;             // 3Fx4 projective camera estimations
  std::vector<double> timings;         // Timings for reconstruction steps
  Eigen::MatrixXd points;              // 4xN projective point estimations
  std::vector<int> pathway;            // Order of views and points addition
  std::vector<std::vector<int>> fixed; // Constraints used for views and points

  Model(int F, int N, int k)
      : inliers(Eigen::MatrixXi::Zero(F, N)),
        cameras(Eigen::MatrixXd::Zero(3 * F, 4)),
        points(Eigen::MatrixXd::Zero(4, N)) {
    // Initialize vectors with default sizes/values as needed
    timings.resize(4, 0.0); // Assuming there are 4 timing steps
    pathway.resize(k, 0);   // Adjust 'k' based on your process
    fixed.resize(F + N);    // Adjust based on your F and N
  }
};

struct SfMModelSeries {
  std::vector<Model> models;

  // Add other fields from `data` structure as needed
  // DataStructure data; // Assuming you have a DataStructure defined elsewhere

  SfMModelSeries(int numModels, int F, int N, int k) {
    models.reserve(numModels);
    for (int i = 0; i < numModels; ++i) {
      models.emplace_back(F, N, k);
    }
  }
};

} // namespace DataStructures
