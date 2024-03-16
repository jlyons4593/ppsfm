#include <Eigen/Dense>
#include <vector>
#pragma once
namespace Options {

inline Eigen::MatrixXi eligibility_view;

const std::vector<int> ELIGIBILITY_POINTS = {10, 9, 8, 7, 6, 5, 4, 3, 2};
const int MAX_LEVEL_POINTS = 7;
const int SCORE_LEVEL = 6;
const int MIN_COMMON_INIT = 200;
const int MAX_MODELS = 1;

const int INIT_LEVEL_POINTS = 5;
const int MAX_LEVEL_VIEWS = 4;
const int INIT_LEVEL_VIEWS = 0;
const bool ROBUST_ESTIMATION = 1;
const double CONFIDENCE = 99.99;
const double RANK_TOLERANCE = 1e-5;
const double SYSTEM_THRESHOLD = 1e-1;
const double OUTLIER_THRESHOLD = 2*1.96;
const int MAX_ITERATION_REFINEMENT=50;
const double MIN_CHANGE_LOCAL_REFINEMENT = 5e-4;
const int MAX_ITER_FINAL_REFINE = 100;
const double MIN_CHANGE_FINAL_REFINE = 1e-4;
const double MIN_CHANGE_GLOBAL_REFINE = 3e-4;
const Eigen::Vector2i MAX_ITERATION_ROBUST(500,1000);
const Eigen::Vector<int, 9> MINIMAL_VIEW{10, 10, 9, 9, 9, 8, 8, 8, 7};
inline Eigen::MatrixXi get_eligibility_view() {
  if (eligibility_view.size() > 0) {
    return eligibility_view;
  } else {
    int data[2][9] = {{84, 60, 48, 36, 24, 18, 12, 10, 8},
                      {8, 4, 4, 2, 2, 2, 1, 1, 1}};
    Eigen::Matrix<int, 2, 9> eligibility_view;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 9; ++j) {
        eligibility_view(i, j) = data[i][j];
      }
    }
    return eligibility_view;
  }
}
}; // namespace Options
