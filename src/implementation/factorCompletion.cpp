#include "factorCompletion.h"
#include "estimate_robust_views.h"
#include "options.h"
#include "pyramidalVisibilityScore.h"
#include "logger.h"
#include <iostream>
#include <numeric>
#include <vector>

void FactorCompletion::process() {

  Logger::logSection("Factor Completion");

  int num_points = data.visible.cols();
  int num_views = data.visible.rows();
  int number_of_visible = (data.visible.array() != 0).count();

  inliers = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero(
      data.visible.rows(), data.visible.cols());
  check_expand_init(num_points, num_views);

  int init_refine = 1;
  int new_init_refine = 1;

  Eigen::Vector2i prev_last_path;
  prev_last_path << 1, 1;

  int last_dir_change = 1;

  rejected_views.resize(num_views);
  rejected_points.resize(num_points);
  rejected_views = Eigen::VectorXi::Zero(num_views);
  rejected_points = Eigen::VectorXi::Zero(num_points);

  init_pvs(num_views);

  bool level_changed = true;
  int level_views =
      std::max(Options::INIT_LEVEL_VIEWS, Options::MAX_LEVEL_VIEWS);
  int level_points =
      std::max(Options::INIT_LEVEL_POINTS, Options::MAX_LEVEL_POINTS);

  int number_of_known_views = camera_variables.negative_pathway.size();
  int number_of_known_points = camera_variables.positive_pathway.size();
  int number_of_added_points = number_of_known_points;
  int number_of_added_views = number_of_known_views;

  int num_iter = 0;
  int iter_refined = 0;

  Logger::logSubsection("Adding Points And Views");

  while (
      (number_of_known_points < num_points ||
       number_of_known_views < num_views) &&
      (level_changed || number_of_added_views + number_of_added_points > 0)) {
    level_changed = false;
    number_of_added_views = 0;
    number_of_added_points = 0;

    if (number_of_known_views < num_views) {
      Eigen::MatrixXi eligibility_view = Options::get_eligibility_view();
      Eigen::VectorXi threshold = eligibility_view.col(level_views);
      std::pair<Eigen::RowVectorXd, std::vector<int>> eligible_scores =  search_eligible_views(threshold, rejected_views);
      Eigen::RowVectorXd scores = eligible_scores.first;
      std::vector<int> eligible_views = eligible_scores.second;
      
      
      if(eligible_views.size()!=0){
        int number_added_views = try_adding_views(eligible_views, level_views);

      }
    }
  }
}

int FactorCompletion::try_adding_views(std::vector<int> eligibles, int level_views){
  int nums_added = 0;
  auto& known_points = camera_variables.positive_pathway;
  // for (int idx =0; idx<eligibles.size(); ++idx){
  for (int idx =0; idx<1; ++idx){

    Eigen::RowVectorXi visible_points(known_points.size());
    for (size_t i = 0; i < known_points.size(); ++i) {
      if (known_points[i] < data.visible.cols()) {
        visible_points(i) = data.visible(eligibles[idx], known_points[i]);
      } else {
        visible_points(i) = 0; 
      }
    }

    Eigen::VectorXi correct_visible_points(known_points.size());
    int count =0;
    for(int i = 0; i < visible_points.size(); ++i) {
      if(visible_points(i)) {
        correct_visible_points(count) = known_points(count);
        count++;
      } 
    }
    if(Options::ROBUST_ESTIMATION){
      EstimatedRobustViews* rob = new EstimatedRobustViews(data, camera_variables,correct_visible_points, eligibles[idx], rejected_views(eligibles[idx]), level_views);
    }else{

    }
    
  }
  return 0;
}
std::pair<Eigen::RowVectorXd ,std::vector<int>> FactorCompletion::search_eligible_views(Eigen::VectorXi thresholds,
                                             Eigen::VectorXi rejected_views) {
  int number_of_views = data.visible.rows();
  // std::cout<<number_of_views<<std::endl;
  Eigen::Array<bool, Eigen::Dynamic, 1> unknown_views =
      Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(number_of_views, 1, true);
  for (int index : camera_variables.negative_pathway) {
    unknown_views(index) = false;
  }
  Eigen::RowVectorXd scores = Eigen::RowVectorXd::Zero(number_of_views);
  Eigen::RowVectorXd num_visible_points =
      Eigen::RowVectorXd::Zero(number_of_views);
  num_visible_points.resize(data.visible.rows());

  for(int i = 0; i < data.visible.rows(); ++i) {
    if(unknown_views[i]) { // Check if the view is unknown
      int sum = 0;
      for(int kp : camera_variables.positive_pathway) {
        if(kp < data.visible.cols()) { // Ensure kp index is within bounds
          sum += data.visible(i, kp);
        }
      }
      num_visible_points(i) = sum; // Assign the sum for this view
    }
  }
  std::vector<bool> eligible_unknown;
  // setting eligible unknown
  for (int i = 0; i < num_visible_points.size(); ++i) {
    if (unknown_views[i]) {
      bool condition1 = num_visible_points[i] > rejected_views[i];
      bool condition2 = num_visible_points[i] >= thresholds[0];
      eligible_unknown.push_back(condition1 && condition2);
    }
  }
  // std::cout<<eligible_unknown[58]<<std::endl;
  for (int i = 0; i < unknown_views.size(); ++i) {
    if (unknown_views[i]) {
      unknown_views[i] = eligible_unknown[i];
    }
  }
  for (size_t i = 0; i < unknown_views.size(); ++i) {
    if (unknown_views(i)) {
      scores(i) = pvs_scores[i]->computeScore(true) * 100;
    }
  }

  std::vector<int> eligibles;
  if ((scores.array() != 0.0).count() > thresholds(1)) {
    eligibles.resize(scores.size());
    std::iota(eligibles.begin(), eligibles.end(),
              0); 
    std::sort(eligibles.begin(), eligibles.end(),
              [&scores](int i1, int i2) { return scores[i1] > scores[i2]; });

    eligibles.resize(thresholds(1));

  }
  else{
    for(int i = 0; i < scores.size(); ++i) {
      if(scores(i) != 0) {
        eligibles.push_back(i); 
      }
    }
  }
  return {scores, eligibles};

}

void FactorCompletion::init_pvs(int num_views) {

  for (int view = 0; view < num_views;  ++view) {
    int rowStart = view * 3;
    int rowEnd = view * 3 + 1;

    std::vector<int> visible_indices;

    for (int kp = 0; kp < camera_variables.positive_pathway.size(); ++kp) {
      if (data.visible(view, camera_variables.positive_pathway[kp])) {
        visible_indices.push_back(camera_variables.positive_pathway[kp]);
      }
    }
    Eigen::MatrixXd projs(2, visible_indices.size()); 

    for(size_t i = 0; i < visible_indices.size(); ++i) {
      projs.col(i) = data.image_measurements.block(rowStart, visible_indices[i], 2, 1); 
    }
    PyramidalVisibilityScore *pvs_score = new PyramidalVisibilityScore(
        image_sizes(0, 0), image_sizes(1, 0), Options::SCORE_LEVEL, projs);

    // If statement needed if multiple models are being calculated
    pvs_scores.push_back(pvs_score);
  }
}

void FactorCompletion::check_expand_init(int num_points, int num_views) {

  auto &pathwayRef = camera_variables.pathway;
  auto &fixedRef = camera_variables.fixed;
  assert(pathwayRef.rows() == 1 && "Initial pathway should be a row");
  assert(pathwayRef.size() == camera_variables.fixed.size());

  last_path = pathwayRef.size();
  int new_zeros_length = num_views + num_points - last_path;

  pathwayRef.conservativeResize(Eigen::NoChange,
                                pathwayRef.cols() + new_zeros_length);
  pathwayRef.tail(new_zeros_length).setZero();
  fixedRef.conservativeResize(fixedRef.size() + new_zeros_length);

  // Initialize the new elements as empty std::vector<int> to mimic empty cells
  for (int i = fixedRef.size() - new_zeros_length; i < fixedRef.size(); ++i) {
    fixedRef(i) = std::vector<int>{};
  }
  int negative_pathway_count = camera_variables.negative_pathway.size();
  int positive_pathway_count = camera_variables.positive_pathway.size();

  // Perform the assertions
  assert(camera_variables.cameras.cols() == 4 &&
         (camera_variables.cameras.rows() == 3 * num_views ||
          camera_variables.cameras.rows() == 3 * negative_pathway_count) &&
         "Bad size for the initial cameras matrix");

  assert(camera_variables.points.rows() == 4 &&
         (camera_variables.points.cols() == num_points ||
          camera_variables.points.cols() == positive_pathway_count) &&
         "Bad size for the initial points matrix");

  if (camera_variables.cameras.rows() != 3 * num_views) {
    Eigen::MatrixXd old_camera = camera_variables.cameras;
    // Fill the matrix with NaN values
    camera_variables.cameras.resize(3 * num_views, 4);
    camera_variables.cameras.fill(std::numeric_limits<double>::quiet_NaN());

    std::vector<double> selectedScaled;
    for (int i = 0; i < negative_pathway_count; ++i) {
      selectedScaled.push_back(3 * camera_variables.negative_pathway[i]);
    }

    Eigen::VectorXi idx(selectedScaled.size() * 3);
    for (int i = 0; i < selectedScaled.size(); i++) {
      for (int k = 0; k < 3; ++k) {
        idx(i * 3 + k) = selectedScaled[i];
      }
    }
    for (int i = 1; i < idx.size(); i += 3) {
      idx[i] += 1;
      if (i + 1 < idx.size()) {
        idx[i + 1] += 2;
      }
    }

    for (int i = 0; i < idx.size(); ++i) {
      if (idx[i] < camera_variables.cameras.rows()) {
        camera_variables.cameras.row(idx[i]) = old_camera.row(i);
      }
    }
  }
  if (camera_variables.points.cols() != num_points) {
    Eigen::MatrixXd old_points = camera_variables.points;

    camera_variables.points.resize(4, num_points);
    camera_variables.points.fill(std::numeric_limits<double>::quiet_NaN());

    for (int i = 0; i < camera_variables.positive_pathway.size(); i++) {

      camera_variables.points.col(camera_variables.positive_pathway(i)) =
          old_points.col(i);
    }
  }

  // Pathway variables different cus -0 is still 0 so using negatives with 0
  // indexing doesnt work
  int number_of_inliers = (inliers.array() > 0).count();
  int number_of_pathway = (camera_variables.pathway.array() > 0).count();

  if (number_of_inliers == 0 && negative_pathway_count == 2) {
    for (int i = 0; i < negative_pathway_count; ++i) {
      for (int j = 0; j < positive_pathway_count; ++j) {
        inliers(camera_variables.negative_pathway(i),
                camera_variables.positive_pathway(j)) = 1;
      }
    }
  }
}
