
#include "dataStructures.hpp"
#include "logger.h"
#include "options.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <numeric>
#include <ostream>
#include <variant>

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
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> visible_pair= filterVisibleMatrix(dense_measurements);

    data.visible = visible_pair.first; 
    dense_measurements = visible_pair.second;
    // Eigen::RowVectorXi visibilitySum = data.visible.cast<int>().colwise().sum();
    // data.removed_points=
    //     visibilitySum.array() <
    //     Options::ELIGIBILITY_POINTS[Options::MAX_LEVEL_POINTS];
    //
    // Logger::logSubsection("Removing Less data.visible Points");
    // std::vector<bool> pointsToRemove(data.visible.cols(), false);
    // for (int i = 0; i < data.removed_points.size(); ++i) {
    //   pointsToRemove[i] = data.removed_points(i);
    // }
    // Eigen::Array<bool, 1, Eigen::Dynamic> keep_pts =
    //     data.removed_points.unaryExpr([](bool v) { return !v; });
    //
    // // Create a list of indices to keep
    // std::vector<int> indices_to_keep;
    // for (int i = 0; i < keep_pts.size(); ++i) {
    //   if (keep_pts(i)) {
    //     indices_to_keep.push_back(i);
    //   }
    // }
    //
    // Eigen::MatrixXd filtered_visible(data.visible.rows(), indices_to_keep.size());
    // Eigen::MatrixXd filtered_orig_meas(dense_measurements.rows(),
    //                                    indices_to_keep.size());
    //
    // // Step 2: Filter columns
    // for (size_t i = 0; i < indices_to_keep.size(); ++i) {
    //   filtered_visible.col(i) = data.visible.col(indices_to_keep[i]);
    //   filtered_orig_meas.col(i) = dense_measurements.col(indices_to_keep[i]);
    // }

    // data.visible = filtered_visible;

    number_of_visible = (data.visible.array() != 0).count();
    number_of_views = data.visible.rows();
    number_of_points = data.visible.cols();
  }

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
  eliminate_dlt(const Eigen::MatrixXd &measurements) {
    Eigen::MatrixXd cpm_meas = measurements.transpose() / measurements.array().square().sum();
    Eigen::MatrixXd data = cpm(measurements);
    Eigen::MatrixXd sub_data = data.topRows(2);
    return {sub_data, cpm_meas};
  }

Eigen::MatrixXd cpm(const Eigen::MatrixXd& v) {
    assert(v.rows() == 3 && v.cols() == 1);

    Eigen::MatrixXd M(3, 3);
    M << 0, -v(2), v(1),
         v(2), 0, -v(0),
        -v(1), v(0), 0;
    return M;
}
  DataStructures::SfMData getData(){
    return data;
  }
  // THIS FUNCTION NOW LOOKS CORRECT
  void process(Eigen::SparseMatrix<double>& measurements) {
    Logger::logSection("Prepare Data");
    dense_measurements = Eigen::MatrixXd(measurements);

    Logger::logSubsection("Creating visibility matrix");

    handleVisibility();
    // ALL VARS CORRECT

    data.image_measurements =
        Eigen::MatrixXd::Zero(3 * number_of_views, number_of_points);
    data.normalised_measurements =
        Eigen::MatrixXd::Zero(3 * number_of_views, number_of_points);
    data.normalisations = Eigen::MatrixXd::Constant(
        3 * number_of_views, 3, std::numeric_limits<double>::quiet_NaN());
    for (int j = 0; j < number_of_views; ++j) {
      Eigen::RowVectorXd visible_points = data.visible.row(j);

      std::vector<int> vis_pts_indices; 
      // Iterate over 'visible_points' to find and store indices of visible
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
  
    Eigen::MatrixXd meas(1,3);
    for (size_t k = 0; k < number_of_visible; ++k) { 
      int startRow = 3 * view_idx[k]; 

      int colIdx =
          point_idx[k]; 
      if (startRow >= 0 && colIdx >= 0 &&
          startRow + 2 < this->data.normalised_measurements.rows() &&
          colIdx < this->data.normalised_measurements.cols()) {
        meas = this->data.normalised_measurements.block(startRow, colIdx, 3, 1);
      }

      auto [cm, pm] = eliminate_dlt(meas);

      pinv_meas_i.segment(3 * k, 3) = Eigen::VectorXd::Constant(3, view_idx[k]);

      int startIdx = 3 * (point_idx[k]); 

      if (startIdx + 2 < pinv_meas_j.size()) {
        for (int i = 0; i < 3; ++i) {
          pinv_meas_j(k*3+i) =
              startIdx + i; 
        }
      }

      for (size_t i = 0; i < pm.size(); ++i) {
        if (3*k + i < pinv_meas_v.size()) {
          pinv_meas_v(k*3 + i) = pm(i);
        }
      }
      int idx1 = 2 * view_idx[k];
      int idx2 = 2 * view_idx[k] + 1;

      // Calculate start index for data_i, data_j, data_v updates

      // Update data_i with idx
      data_i.segment<2>(k*6, 2) << idx1, idx2;
      data_i.segment<2>(k*6 + 2, 2) << idx1, idx2;
      data_i.segment<2>(k*6 + 4, 2) << idx1, idx2;

      data_j(k*6) = 3 * point_idx[k];     
      data_j(k*6 + 1) = 3 * point_idx[k]; 
      data_j(k*6 + 2) = 3 * point_idx[k] + 1;
      data_j(k*6 + 3) = 3 * point_idx[k] + 1;
      data_j(k*6 + 4) = 3 * point_idx[k] + 2;
      data_j(k*6 + 5) = 3 * point_idx[k] + 2;

      std::vector<double> cm_flattened;
      for (int j = 0; j < cm.cols(); ++j) {
        for (int i = 0; i < cm.rows(); ++i) {
          cm_flattened.push_back(cm(i, j));
        }
      }
      // data_v
      for (size_t i = 0; i < cm_flattened.size(); ++i) {
        if (k*6 + i < data_v.size()) {
          data_v(k*6 + i) = cm_flattened[i];
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

     rows = static_cast<int>(pinv_meas_i.maxCoeff()) + 1;
     cols = static_cast<int>(pinv_meas_j.maxCoeff()) + 1;
    
     // Initialize the data matrix with zeros
     data.pseudo_inverse_measurements = Eigen::MatrixXd::Zero(rows, cols);
    
     for(int p = 0; p < pinv_meas_i.size(); ++p) {
         int row = static_cast<int>(pinv_meas_i(p));
         int col = static_cast<int>(pinv_meas_j(p));
         data.pseudo_inverse_measurements(row, col) = pinv_meas_v(p);
     }

  }

std::pair<Eigen::MatrixXd,Eigen::MatrixXd> filterVisibleMatrix(Eigen::MatrixXd &dense_measurements) {

    // creating a matrix of same size as dense measurements that is a
    // binary visibility matrix
    Eigen::MatrixXd visible = (dense_measurements.array() != 0).cast<double>();
    

    // filter rows to ensure visibility in both subsequent rows
    // Eigen::MatrixXd filtered_visible(visible.rows() / 2, visible.cols());
    // for (int i = 0; i < visible.rows() / 2; ++i) {
    //   filtered_visible.row(i) =
    //       visible.row(2 * i).cwiseProduct(visible.row(2 * i + 1));
    // }
    // visible = filtered_visible;

    Eigen::MatrixXd reduced(visible.rows() / 2, visible.cols());
    for (int i = 0; i < visible.rows(); i += 2) {
        for (int j = 0; j < visible.cols(); ++j) {
            reduced(i / 2, j) = visible(i, j) && visible(i + 1, j);
        }
    }

    visible = reduced;
    // Compute the sum of each column to count the number of views each point is
    // data.visible in

    Eigen::VectorXd colSums = visible.colwise().sum();

    // Determine which columns have a sum less than the specified threshold
    Eigen::Array<bool, 1, Eigen::Dynamic> rm_pts = colSums.array() < Options::ELIGIBILITY_POINTS[Options::MAX_LEVEL_POINTS];
    Eigen::Array<bool, 1, Eigen::Dynamic> keep_pts = !rm_pts;

    int numColsToKeep = keep_pts.count();

    // Create a new matrix to hold the filtered columns
    Eigen::MatrixXd filtered_visible(visible.rows(), numColsToKeep);
    Eigen::MatrixXd filtered_meas(dense_measurements.rows(), numColsToKeep);

    // Iterate over keep_pts and copy the columns to keep into the new matrix
    int colIndex = 0;
    for (int i = 0; i < keep_pts.size(); ++i) {
        if (keep_pts(i)) {
            filtered_visible.col(colIndex) = visible.col(i);
            filtered_meas.col(colIndex) = dense_measurements.col(i);
            colIndex++;
        }
    }

    // Update data.visible with the filtered matrix
    return {filtered_visible,filtered_meas};
  }
  //
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
  normtrans(const Eigen::MatrixXd &inputPoints, bool isotropic = true) {
    Eigen::MatrixXd points = inputPoints;
    int number_of_points = points.cols();
    int dimension = points.rows();

    bool homogeneous = (points.row(dimension - 1).array() == 1).all();
    if (homogeneous) {
      --dimension; 
                  
      points.conservativeResize(
          dimension,
          Eigen::NoChange); 
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
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
  eliminate_pinv(const Eigen::MatrixXd &measurements) {
    // sparse matrices Look into changing this to keeping this sparse for
    // optimisations later

    // Transpose the measurements and divide by the sum of its squared elements
    Eigen::MatrixXd pseudo_inverse_measurements =
        measurements.transpose() / measurements.array().square().sum();

    Eigen::MatrixXd data = measurements * pseudo_inverse_measurements;

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
