#include "ppsfm_initialisation.h"
#include "logger.h"
#include <iostream>

void Initialisation::process() {

  Logger::logSection("Starting PPSFM Initialisation");

  Logger::logSubsection("Creating Matrices for Computation");

  Eigen::MatrixXi estimated_pairs =
      Eigen::MatrixXi::Zero(view_pairs.rows(), view_pairs.cols());

  Eigen::MatrixXi not_estimated_pairs =
      (estimated_pairs.array() == 0).cast<int>();
  Eigen::VectorXi both_unestimated(not_estimated_pairs.rows());

  Eigen::MatrixXi visibleLogical = visible.cast<int>();

  for (int i = 0; i < not_estimated_pairs.rows(); i++) {
    if (not_estimated_pairs(i, 0) && not_estimated_pairs(i, 1)) {
      both_unestimated(i) = 1;
    }
  }

  Eigen::VectorXi one_estimated(estimated_pairs.rows());

  for (int i = 0; i < estimated_pairs.rows(); i++) {
    bool any_estimated =
        false; 
    bool all_estimated =
        true; 

    for (int j = 0; j < estimated_pairs.cols(); j++) {
      if (estimated_pairs(i, j) == 0) {
        all_estimated = false; 
      } else {
        any_estimated =
            true; 
      }
    }
    one_estimated(i) = any_estimated && !all_estimated;
  }

  std::vector<int> sorted_idx;

  // Find indices of non-zero elements in both_unestimated
  for (int i = 0; i < both_unestimated.size(); i++) {
    if (both_unestimated(i) != 0) {
      sorted_idx.push_back(i);
    }
  }

  // Find indices of non-zero elements in one_estimated and append
  for (int i = 0; i < one_estimated.size(); i++) {
    if (one_estimated(i) != 0) {
      sorted_idx.push_back(i);
    }
  }

  // MAIN WORK LOOP

  Logger::logSubsection("Main PPSFM Initialisation Loop");

  for (int i = 0; i < sorted_idx.size(); i++) {
    // for(int i=0; i<1; i++){
    int first_view =
        view_pairs(sorted_idx[i], 0); 
    int second_view = view_pairs(sorted_idx[i], 1);

    // logical and for creating visible points rowVector
    auto logicalAnd = [](int a, int b) -> int { return a && b; };
    Eigen::VectorXi visible_points =
        visible.row(first_view)
            .binaryExpr(visible.row(second_view), logicalAnd);
    // std::cout<<visible_points.size()<<std::endl;
    std::vector<int> first_view_idx;
    std::vector<int> second_view_idx;
    for (int i = 3 * first_view; i < 3 * first_view + 2; i++) {
      first_view_idx.push_back(i);
    }
    for (int i = 3 * second_view; i < 3 * second_view + 2; i++) {
      second_view_idx.push_back(i);
    }
    int number_of_visible = (visible_points.array() != 0).count();

    Eigen::MatrixXd block1(second_view_idx.size(), number_of_visible);
    Eigen::MatrixXd block2(first_view_idx.size(), number_of_visible);

    // Fill measurement block 1
    for (int i = 0; i < second_view_idx.size(); ++i) {
      for (int j = 0; j < visible_points.size(); ++j) {
        if (visible_points(j) > 0) {
          block1(i, j) = measurement(second_view_idx[i], j);
        }
      }
    }
    // Fill measuremnent block 2
    for (int i = 0; i < first_view_idx.size(); ++i) {
      for (int j = 0; j < visible_points.size(); ++j) {
        if (visible_points(j) > 0) {
          block2(i, j) = measurement(first_view_idx[i], j);
        }
      }
    }

    // Call to estimated fundamental matrix
    EstimatedFundamentalMatrix estimated_fundamental_matrix(block1, block2,
                                                            99.99, 1000, 1e-3);

    Eigen::Array<bool, Eigen::Dynamic, 1> inliers =
        estimated_fundamental_matrix.getInliers();
    Eigen::MatrixXd fundamental_matrix =
        estimated_fundamental_matrix.getFundamentalMatrix();

    int status = 0;

    if (status == 0) {
      std::vector<int> initial_views = {first_view, second_view};

      std::vector<int> visible_idx;
      for (int i = 0; i < visible_points.size(); ++i) {
        if (visible_points(i)) {
          visible_idx.push_back(
              i); // Store zero-based indices of visible points
        }
      }

      std::vector<int> initial_points;
      for (int i = 0; i < inliers.size(); ++i) {
        if (inliers(i)) { // If the inlier_idx at position i is true, select the
                          // corresponding index
          initial_points.push_back(visible_idx[i]);
        }
      }
      camera_variables =
          compute_cams(initial_views, initial_points, fundamental_matrix);
      break;
    }
  }
  Logger::logSection("End of Camera Var Computation");
}

DataStructures::ComputedCameraPoints
Initialisation::compute_cams(std::vector<int> initial_views,
                             std::vector<int> initial_points,
                             Eigen::MatrixXd fundamental_matrix) {
  DataStructures::ComputedCameraPoints camera_points;
  Eigen::Vector3i first_view(3 * initial_views[0], 3 * initial_views[0] + 1,
                             3 * initial_views[0] + 2);
  Eigen::Vector3i second_view(3 * initial_views[1], 3 * initial_views[1] + 1,
                              3 * initial_views[1] + 2);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(fundamental_matrix,
                                        Eigen::ComputeFullV);
  Eigen::MatrixXd V = svd.matrixV();

  // Epipole is the 3rd column
  Eigen::Vector3d epipole = V.col(2);

  Eigen::MatrixXd proj_depths = Eigen::MatrixXd::Ones(2, initial_points.size());

  for (int point_idx = 0; point_idx < initial_points.size(); point_idx++) {
    Eigen::Vector3d vector_from_measurement;
    for (int i = 0; i < 3; ++i) {
      vector_from_measurement(i) =
          measurement(second_view[i], initial_points[point_idx]);
    }

    // Compute the cross product
    Eigen::VectorXd point_epipole = vector_from_measurement.cross(epipole);
    int col_idx = initial_points[point_idx];
    Eigen::Vector3d meas_vec;
    for (int i = 0; i < 3; ++i) {
      meas_vec(i) = measurement(first_view(i), col_idx);
    }

    double numerator = (meas_vec.transpose() * fundamental_matrix *
                        point_epipole)(0, 0); // This extracts the scalar value
    double denominator =
        (point_epipole.transpose() * point_epipole)(0, 0); // Also a scalar

    double value = std::abs(numerator / denominator);
    proj_depths(1, point_idx) = value;
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
  camera_points.positive_pathway.resize(initial_points.size());
  camera_points.negative_pathway.resize(initial_views.size());
  camera_points.positive_pathway(0) = initial_points[0];
  camera_points.negative_pathway(0) = initial_views[0];
  pathway.push_back(initial_points[0]);
  pathway.push_back(-initial_views[0]);


  camera_points.positive_pathway(1) = initial_points[1];
  camera_points.negative_pathway(1) = initial_views[1];
  pathway.push_back(initial_points[1]);
  pathway.push_back(-initial_views[1]);


  for (size_t i = 2; i < initial_points.size(); ++i) {
    pathway.push_back(initial_points[i]);
    camera_points.positive_pathway(i) = initial_points[i];
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
  for (auto i : first_view) {
    combined_indices.push_back(i);
  }
  combined_indices.insert(combined_indices.end(), second_view.begin(),
                          second_view.end());

  // Create a scaled_measurement matrix with the correct dimensions
  Eigen::MatrixXd scaled_measurement(combined_indices.size(),
                                     initial_points.size());

  // Fill in scaled_measurement by manually selecting the appropriate elements
  for (size_t i = 0; i < combined_indices.size(); ++i) {
    for (size_t j = 0; j < initial_points.size(); ++j) {
      scaled_measurement(i, j) =
          measurement(combined_indices[i], initial_points[j]);
    }
  }
  assert(proj_depths.rows() == 2);
  //
  Eigen::MatrixXd adjusted_proj_depths(6, proj_depths.cols());
  for (int i = 0; i < proj_depths.cols(); ++i) {
    for (int j = 0; j < 3;
         ++j) { // First 3 rows are copies of the first row of proj_depths
      adjusted_proj_depths.row(j)(i) = proj_depths.row(0)(i);
    }
    for (int j = 3; j < 6;
         ++j) { // Next 3 rows are copies of the second row of proj_depths
      adjusted_proj_depths.row(j)(i) = proj_depths.row(1)(i);
    }
  }
  //
  // Perform element-wise multiplication
  scaled_measurement =
      scaled_measurement.array() * adjusted_proj_depths.array();
  Eigen::JacobiSVD<Eigen::MatrixXd> scaled_svd(
      scaled_measurement, Eigen::ComputeFullV | Eigen::ComputeThinU);

  // Extract the matrices U, S, and V
  Eigen::MatrixXd scaled_U = scaled_svd.matrixU();
  Eigen::MatrixXd scaled_S = scaled_svd.singularValues().asDiagonal();
  Eigen::MatrixXd scaled_V = scaled_svd.matrixV();

  Eigen::VectorXd sqrt_singularValues =
      scaled_svd.singularValues().array().sqrt();

  Eigen::MatrixXd sqrt_S_4 = sqrt_singularValues.head(4).asDiagonal();

  // Now perform the multiplication with the first four columns of U
  camera_points.cameras = scaled_U.leftCols(4) * sqrt_S_4;
  camera_points.points = (sqrt_S_4 * scaled_V.leftCols(4).transpose());

  // HACKY CONVERSION PLEASE FIX
  Eigen::Map<Eigen::VectorXi> temp_vector(pathway.data(), pathway.size());
  camera_points.pathway = temp_vector;
  return camera_points;
}

DataStructures::ComputedCameraPoints Initialisation::getCameraPoints() {
  return camera_variables;
}

Initialisation::Initialisation(const Eigen::MatrixXd &measurement,
                               const Eigen::MatrixXd &visible,
                               const Eigen::MatrixXd &view_pairs,
                               const Eigen::VectorXd &affinity)
    : measurement(measurement), visible(visible), view_pairs(view_pairs),
      affinity(affinity) {}
