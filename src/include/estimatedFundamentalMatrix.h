#include "dataStructures.hpp"
#include <Eigen/Dense>
#include <algorithm> // For std::shuffle
#include <iostream>
#include <numeric> // For std::iota
#include <ostream>
#include <random> // For std::default_random_engine
#include <vector>

class EstimatedFundamentalMatrix {
private:
  Eigen::MatrixXd fundamental_matrix;
  Eigen::Array<bool, Eigen::Dynamic, 1> inliers;
  std::random_device rd;
  std::default_random_engine rng;

  Eigen::VectorXd linvec(const Eigen::VectorXd &a, const Eigen::VectorXd &b);

  std::vector<int> get_random_sequence(int number_of_projections);

  Eigen::MatrixXd estimate(Eigen::MatrixXd coeffs, bool enforce_rank);

  std::pair<Eigen::MatrixXd, Eigen::Array<bool, Eigen::Dynamic, 1>>
  estimateRansac(Eigen::MatrixXd coeffs, double confidence, int max_iter,
                 double dist_thresh);

  Eigen::MatrixXd applyLinvecToProjections(const Eigen::MatrixXd &projs1,
                                           const Eigen::MatrixXd &projs2);

public:
  Eigen::MatrixXd getFundamentalMatrix() { return fundamental_matrix; }
  Eigen::Array<bool, Eigen::Dynamic, 1> getInliers() { return inliers; }
  EstimatedFundamentalMatrix(
      const Eigen::MatrixXd &projs1, // Projections in the first image (2xN)
      const Eigen::MatrixXd
          &projs2,       // Corresponding projections in the second image (2xN)
      double confidence, // Confidence to stop robust estimation early in RANSAC
      int max_iter,      // Maximum number of iterations in RANSAC
      double dist_thresh // Distance threshold for computing inliers in RANSAC
      )
      : rng(rd()) {

    // Setting default random device
    //
    Eigen::MatrixXd coeffs = applyLinvecToProjections(projs1, projs2);
    // Coeffs are correct

    // Currently only works for Ransac Method
    std::pair<Eigen::MatrixXd, Eigen::Array<bool, Eigen::Dynamic, 1>>
        fund_mat_inliers =
            estimateRansac(coeffs, confidence, max_iter, dist_thresh);

    fundamental_matrix = fund_mat_inliers.first;
    inliers = fund_mat_inliers.second;
  }
};
