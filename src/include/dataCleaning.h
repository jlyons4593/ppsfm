
#include "dataStructures.hpp"
#include "logger.h"
#include "options.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <numeric>
#include <ostream>

#pragma once
class DataCleaningStage {
private:
  DataStructures::SfMData data;
  Eigen::MatrixXd dense_measurements;
  int number_of_views;
  int number_of_points;
  int number_of_visible;
  Eigen::VectorXd data_i;
  Eigen::VectorXd data_j;
  Eigen::VectorXd data_v;
  Eigen::VectorXd pinv_meas_i;
  Eigen::VectorXd pinv_meas_j;
  Eigen::VectorXd pinv_meas_v;


  void sizeDataPinvVectors(){

    this->data_i.resize(6*number_of_visible);
    this->data_j.resize(6*number_of_visible);
    this->data_v.resize(6*number_of_visible);
    this->pinv_meas_i.resize(3*number_of_visible);
    this->pinv_meas_j.resize(3*number_of_visible);
    this->pinv_meas_v.resize(3*number_of_visible);

    this->data_i.setZero();
    this->data_j.setZero();
    this->data_v.setZero();
    this->pinv_meas_i.setZero();
    this->pinv_meas_j.setZero();
    this->pinv_meas_v.setZero();
  }

public:
   void handleVisibility(){


    Logger::logSubsection("Creating visibility matrix");
    data.visible = filterVisibleMatrix(dense_measurements);
    Eigen::RowVectorXi visibilitySum = data.visible.cast<int>().colwise().sum();
    data.removed_points=
        visibilitySum.array() <
        Options::ELIGIBILITY_POINTS[Options::MAX_LEVEL_POINTS];

    Logger::logSubsection("Removing Less data.visible Points");
    std::vector<bool> pointsToRemove(data.visible.cols(), false);
    for (int i = 0; i < data.removed_points.size(); ++i) {
      pointsToRemove[i] = data.removed_points(i);
    }
    Eigen::Array<bool, 1, Eigen::Dynamic> keep_pts =
        data.removed_points.unaryExpr([](bool v) { return !v; });

    // Create a list of indices to keep
    std::vector<int> indices_to_keep;
    for (int i = 0; i < keep_pts.size(); ++i) {
      if (keep_pts(i)) {
        indices_to_keep.push_back(i);
      }
    }

    // Assuming 'data.visible' and 'orig_meas' are Eigen matrices
    Eigen::MatrixXd filtered_visible(data.visible.rows(), indices_to_keep.size());
    Eigen::MatrixXd filtered_orig_meas(dense_measurements.rows(),
                                       indices_to_keep.size());

    // Step 2: Filter columns
    for (size_t i = 0; i < indices_to_keep.size(); ++i) {
      filtered_visible.col(i) = data.visible.col(indices_to_keep[i]);
      filtered_orig_meas.col(i) = dense_measurements.col(indices_to_keep[i]);
    }
    data.visible = filtered_visible;
    dense_measurements = filtered_orig_meas;

    number_of_visible = (data.visible.array() != 0).count();
    number_of_views = data.visible.rows();
    number_of_points = data.visible.cols();
  }

  DataStructures::SfMData getData(){
    return data;
  }
  // THIS FUNCTION NOW LOOKS CORRECT
  void process(Eigen::SparseMatrix<double> &measurements) {
    Logger::logSection("Prepare Data");
    dense_measurements = Eigen::MatrixXd(measurements);

    Logger::logSubsection("Creating visibility matrix");

    handleVisibility();
    // ALL VARS CORRECT

    data.image_measurements =
        Eigen::MatrixXd::Zero(3 * number_of_views, number_of_points);
    // printColsRows(image_measurements, "Image Measurements");
    data.normalised_measurements =
        Eigen::MatrixXd::Zero(3 * number_of_views, number_of_points);
    // printColsRows(normalised_measurements, "Normalised Measurements");
    data.normalisations = Eigen::MatrixXd::Constant(
        3 * number_of_views, 3, std::numeric_limits<double>::quiet_NaN());
    // printColsRows(normalisations, "Normalisations");
    for (int j = 0; j < number_of_views; ++j) {
      Eigen::RowVectorXd visible_points = data.visible.row(j);

      std::vector<int> vis_pts_indices; // This will need to be filled with
      // Iterate over 'visible_points' to find and store indices of visible
      // points
      for (int i = 0; i < visible_points.size(); ++i) {
        if (visible_points(i) > 0) { 
          vis_pts_indices.push_back(i);
        }
      }

      for (int idx : vis_pts_indices) {
        data.image_measurements.block(j * 3, idx, 2, 1) =
            dense_measurements.block((j) * 2, idx, 2, 1);
      }
      for (int idx : vis_pts_indices) {
        data.image_measurements((j * 3) + 2, idx) = 1; 
      }

      Eigen::MatrixXd block(3, vis_pts_indices.size());
      for (size_t colIndex = 0; colIndex < vis_pts_indices.size(); ++colIndex) {
        int visCol = vis_pts_indices[colIndex];
        block.col(colIndex) = data.image_measurements.block(j * 3, visCol, 3, 1);
      }

      // Call normtrans
      auto [transform, transformed_measurements] = normtrans(block);

      data.normalisations.block(j * 3, 0, 3, data.normalisations.cols()) = transform;
      // Update normalisations
      // EXPLAIN WHAT THIS IS DOING
      for (int i = 0; i < vis_pts_indices.size(); ++i) {
        int point_col = vis_pts_indices[i]; 
                                           

        data.normalised_measurements.block(j * 3, point_col, 3, 1) =
            transformed_measurements.col(i);
      }
    }


    sizeDataPinvVectors();

    std::vector<int> view_idx, point_idx;

    for (int i = 0; i < this->data.visible.rows(); ++i) {
      for (int j = 0; j < this->data.visible.cols(); ++j) {
        if (this->data.visible(i, j) > 0) {
          view_idx.push_back(i);
          point_idx.push_back(j);
        }
      }
    }
    // std::cout<<"view idx size: "<<view_idx.size()<<std::endl;
    // std::cout<<"point idx size: "<<point_idx.size()<<std::endl;
    //

    Eigen::MatrixXd meas;
    for (size_t k = 0; k < number_of_visible; ++k) { 
      int startRow = 3 * view_idx[k] - 3; 

      int colIdx =
          point_idx[k]; 
      // Ensure indices are within bounds
      if (startRow >= 0 && colIdx >= 0 &&
          startRow + 2 < this->data.normalised_measurements.rows() &&
          colIdx < this->data.normalised_measurements.cols()) {
        meas = this->data.normalised_measurements.block(startRow, colIdx, 3, 1);
        // Optionally, use 'meas' here, e.g., print or process it
      }

      auto [cm, pm] = eliminate_pinv(meas);

      pinv_meas_i.segment(3 * k, 3) = Eigen::VectorXd::Constant(3, view_idx[k]);
      int startIdx = 3 * (point_idx[k]); 

      // Assuming pinv_meas_j has been sized correctly
      if (startIdx + 2 < pinv_meas_j.size()) {
        for (int i = 0; i < 3; ++i) {
          pinv_meas_j(startIdx + i) =
              startIdx + i; 
        }
      }

      // Assign values from pm to pinv_meas_v
      for (size_t i = 0; i < pm.size(); ++i) {
        if (k + i < pinv_meas_v.size()) {
          pinv_meas_v(k + i) = pm(i);
        }
      }
      int idx1 = 2 * view_idx[k];
      int idx2 = 2 * view_idx[k] + 1;

      // Calculate start index for data_i, data_j, data_v updates

      // Update data_i with idx
      data_i.segment(k, 2) << idx1, idx2;
      data_i.segment(k + 2, 2) << idx1, idx2;
      data_i.segment(k + 4, 2) << idx1, idx2;

      data_j(k) = 3 * point_idx[k];     
      data_j(k + 1) = 3 * point_idx[k]; 
      data_j(k + 2) = 3 * point_idx[k] + 1;
      data_j(k + 3) = 3 * point_idx[k] + 1;
      data_j(k + 4) = 3 * point_idx[k] + 2;
      data_j(k + 5) = 3 * point_idx[k] + 2;

      std::vector<double> cm_flattened;
      for (int i = 0; i < cm.rows(); ++i) {
        for (int j = 0; j < cm.cols(); ++j) {
          cm_flattened.push_back(cm(i, j));
        }
      }
      // Assuming cm_flattened now contains the flattened values, assign them to
      // data_v
      for (size_t i = 0; i < cm_flattened.size(); ++i) {
        if (k + i < data_v.size()) {
          data_v(k + i) = cm_flattened[i];
        }
      }
    }


    int rows = static_cast<int>(data_i.maxCoeff()) + 1;
    int cols = static_cast<int>(data_j.maxCoeff()) + 1;

    // Initialize the data matrix with zeros
     data.cost_function_data = Eigen::MatrixXd::Zero(rows, cols);

     for(int k = 0; k < data_i.size(); ++k) {
         int row = static_cast<int>(data_i(k));
         int col = static_cast<int>(data_j(k));
         data.cost_function_data(row, col) = data_v(k);
     }
//
    DataStructures::printColsRows(data.cost_function_data, "Cost Function Data");
//
     rows = static_cast<int>(pinv_meas_i.maxCoeff()) + 1;
     cols = static_cast<int>(pinv_meas_j.maxCoeff()) + 1;
    
     // Initialize the data matrix with zeros
     data.pseudo_inverse_measurements = Eigen::MatrixXd::Zero(rows, cols);
    
     for(int k = 0; k < pinv_meas_i.size(); ++k) {
         int row = static_cast<int>(pinv_meas_i(k));
         int col = static_cast<int>(pinv_meas_j(k));
         data.pseudo_inverse_measurements(row, col) = pinv_meas_v(k);
     }
    DataStructures::printColsRows(data.pseudo_inverse_measurements,"pinv_meas");
  }

  Eigen::MatrixXd filterVisibleMatrix(Eigen::MatrixXd &dense_measurements) {

    // Step 1: creating a matrix of same size as dense measurements that is a
    // binary visibility matrix
    Eigen::MatrixXd visible = (dense_measurements.array() != 0).cast<double>();

    // Step 2: Filter rows to ensure visibility in both subsequent rows
    Eigen::MatrixXd filtered_visible(visible.rows() / 2, visible.cols());
    for (int i = 0; i < visible.rows() / 2; ++i) {
      filtered_visible.row(i) =
          visible.row(2 * i).cwiseProduct(visible.row(2 * i + 1));
    }
    visible = filtered_visible;

    // Compute the sum of each column to count the number of views each point is
    // data.visible in
    Eigen::VectorXi visibility_count = visible.cast<int>().colwise().sum();

    int threshold = Options::ELIGIBILITY_POINTS[Options::MAX_LEVEL_POINTS];

    std::vector<int> points_to_keep;
    for (int i = 0; i < visibility_count.size(); ++i) {
      if (visibility_count[i] >= threshold) {
        points_to_keep.push_back(i);
      }
    }
    // Filter data.visible and orig_meas matrices to keep only the columns for points
    // above the threshold
    Eigen::MatrixXd refiltered_visible(visible.rows(), points_to_keep.size());
    for (size_t i = 0; i < points_to_keep.size(); ++i) {
      refiltered_visible.col(i) = visible.col(points_to_keep[i]);
    }

    // Update data.visible with the filtered matrix
    return refiltered_visible;
  }
  //
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
  normtrans(const Eigen::MatrixXd &inputPoints, bool isotropic = true) {
    Eigen::MatrixXd points = inputPoints;
    int number_of_points = points.cols();
    int dimension = points.rows();

    // Check for homogeneity
    bool homogeneous = (points.row(dimension - 1).array() == 1).all();
    if (homogeneous) {
      --dimension; // Adjust for homogeneity by ignoring the last row for
                   // centroid and scaling calculations
      points.conservativeResize(
          dimension,
          Eigen::NoChange); // Temporarily adjust points matrix for processing
    }

    // Calculate centroid
    Eigen::VectorXd centroid = points.rowwise().mean();

    // Calculate difference from centroid
    Eigen::MatrixXd diff = points.colwise() - centroid;

    // Calculate scaling factor
    double scale =
        sqrt(2) / (diff.array().square().colwise().sum().sqrt().mean());

    // Construct transformation matrix
    Eigen::MatrixXd trans = Eigen::MatrixXd::Identity(
        dimension + (homogeneous ? 1 : 0), dimension + (homogeneous ? 1 : 0));
    trans.topLeftCorner(dimension, dimension) *= scale;
    trans.topRightCorner(dimension, 1) = -centroid * scale;
    if (!homogeneous) {
      // If not homogeneous, ensure the last row is appropriate for
      // non-homogeneous coordinates
      trans.conservativeResize(dimension + 1, dimension + 1);
      trans.row(dimension).setZero();
      trans(dimension, dimension) = 1;
    }

    // Transform points
    Eigen::MatrixXd transformedPoints = trans * inputPoints;

    return {trans, transformedPoints};
  }
  /*
   * This will require implementation of the CPM function for signal processing
   * from matlab into C++
   * */
  // std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eliminate_dlt(const
  // Eigen::SparseMatrix<double> &measurements) {
  //   // Calculate cpm measurements
  //   Eigen::SparseMatrix<double> meas_transposed = measurements.transpose();
  //   double sumOfSquares = 0.0;
  //   for (int k = 0; k < measurements.outerSize(); ++k) {
  //       for (Eigen::SparseMatrix<double>::InnerIterator it(measurements, k);
  //       it; ++it) {
 //           sumOfSquares += std::pow(it.value(), 2);
  //       }
  //   }
  //   Eigen::SparseMatrix<double> cpm_meas = meas_transposed / sumOfSquares;
  //
  //   return
  // }

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
  eliminate_pinv(const Eigen::MatrixXd &measurements) {
    // First, convert meas toe dense for operations that are not well-defined on
    // sparse matrices Look into changing this to keeping this sparse for
    // optimisations later

    // Transpose the measurements and divide by the sum of its squared elements
    Eigen::MatrixXd pseudo_inverse_measurements =
        measurements.transpose() / measurements.array().square().sum();

    // Multiply measurements by pseudo_inverse_measurements to get data
    // Since this operation likely results in a dense matrix, we perform it on
    // the dense version of meas
    Eigen::MatrixXd data = measurements * pseudo_inverse_measurements;

    // Adjust specific diagonal elements of data
    if (data.rows() >= 3 && data.cols() >= 3) {
      data(0, 0) -= 1;
      data(1, 1) -= 1;
      data(2, 2) -= 1;
    }

    // Resize data to keep only the first two rows
    Eigen::MatrixXd resizedData(2, data.cols());
    if (data.rows() >= 2) {
      resizedData = data.topRows(2);
    }

    return {resizedData, pseudo_inverse_measurements};
  }
  DataCleaningStage() {}
};
