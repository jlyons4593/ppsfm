#include "factorCompletion.h"
#include "diagnosis.h"
#include "Refinement.h"
#include "estimatedViews.h"
#include "estimatedPoints.h"
#include "estimate_robust_views.h"
#include "estimatedRobustPoints.h"
#include "options.h"
#include "pyramidalVisibilityScore.h"
#include "logger.h"
#include <exception>
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
  // std::cout<<camera_variables.negative_pathway<<std::endl;
  int number_of_known_points = camera_variables.positive_pathway.size();
  int number_of_added_points = number_of_known_points;
  int number_of_added_views = number_of_known_views;

  int num_iter = 0;
  int iter_refine = 0;

  Diagnosis* diag= new Diagnosis(inliers, data, camera_variables, num_iter);
  delete(diag);
  Logger::logSubsection("Adding Points And Views");

  int refinement_count = 0;
  int actual_iter = 0;

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
          // std::cout<<"before search eligible views"<<std::endl;
      std::pair<Eigen::RowVectorXd, std::vector<int>> eligible_scores =  search_eligible_views(threshold, rejected_views);
          // std::cout<<"after search eligible views"<<std::endl;
      Eigen::RowVectorXd scores = eligible_scores.first;
      std::vector<int> eligible_views = eligible_scores.second;
      
      
      if(eligible_views.size()!=0){
        number_of_added_views= try_adding_views(eligible_views, level_views);

          // std::cout<<"after try adding views"<<std::endl;
        if(number_of_added_views>0){
          number_of_known_views = number_of_known_views+number_of_added_views;
          level_points = std::max(1, level_points-1);

          // std::cout<<"before Refinement"<<std::endl;

          refinement_count++;
          Refinement* refinement = new Refinement(data, camera_variables, inliers,init_refine, last_path, false,0 );
          // break;
          // std::cout<<"after Refinement"<<std::endl;
          camera_variables.cameras= refinement->getCameras();
          camera_variables.points= refinement->getPoints();
          delete(refinement);
          Diagnosis* diag = new Diagnosis(inliers, data, camera_variables, num_iter);
          delete diag;


          if(last_dir_change==2){
            prev_last_path[0] = last_path - number_of_added_views;
            last_dir_change =1;
            num_iter++;
          }
        }
      }
      if(!number_of_added_views && level_views<Options::MAX_LEVEL_VIEWS){
        std::cout<<"If statement 1 "<<std::endl;
        level_views++;
        level_changed = true;
      }
    }
    // if(actual_iter==13){
    //   std::cout<<"wrong route"<<std::endl;
    //   std::cout<<camera_variables.negative_pathway.size()<<std::endl;
    // std::cout<<camera_variables.negative_pathway<<std::endl;
    //   break;
    //
    // }
    Eigen::VectorXi eligible_points = search_eligible_points(Options::ELIGIBILITY_POINTS[level_points], rejected_points, actual_iter);

    if(eligible_points.size()!=0){
// std::cout<<camera_variables.cameras.col(0)<<std::endl;
//   throw std::exception();
      std::pair<int, Eigen::VectorXi>added_vars =  try_adding_points(eligible_points, level_views);
        std::cout<<"iut of try adding pts"<<std::endl;
      number_of_added_points = added_vars.first;
      Eigen::VectorXi added= added_vars.second;
        std::cout<<"vars retreaved"<<std::endl;

      if(number_of_added_points>0){
        std::cout<<"into nex if"<<std::endl;
        
        number_of_known_points = number_of_known_points+ number_of_added_points;
        update_pvs(eligible_points(added));

      }
      
      

    }
    if(!number_of_added_points && level_points<Options::MAX_LEVEL_POINTS-Options::DIFFER_LAST_LEVEL){
        std::cout<<"If statement 2 "<<std::endl;
      level_points++;
      level_changed = true;
    }
    if(new_init_refine != std::min(prev_last_path[0],prev_last_path[1])){
        std::cout<<"If statement 3 "<<std::endl;
      init_refine = new_init_refine;
      new_init_refine = std::min(prev_last_path[0], prev_last_path[1]);
    }
    if((num_iter-iter_refine)>Options::GLOBAL_REFINE){
      std::cout<<"If statement 4"<<std::endl;
      Refinement* refinement = new Refinement(data, camera_variables, inliers,init_refine, last_path, false,1 );
      camera_variables.cameras= refinement->getCameras();
      camera_variables.points= refinement->getPoints();
      iter_refine = num_iter;
    }
    actual_iter++;

    std::cout<<"actual iteration number: "<<actual_iter<<std::endl;
    // std::cout<<"max complete iter num: "<<ite<<std::endl;
  } // End of loop
  camera_variables.pathway = camera_variables.pathway.segment(0,last_path);
  std::vector<std::vector<int>> sectionVector;
  sectionVector.assign(camera_variables.fixed.begin(), camera_variables.fixed.begin() + last_path);
  camera_variables.fixed = sectionVector;
  std::cout<<"refinement count: "<<refinement_count<<std::endl;
}
void FactorCompletion::update_pvs(Eigen::VectorXi added_points){
  std::cout<<"updating"<<std::endl;
  // Eigen::MatrixXd visible_added(data.visible.rows(), added_points.size());
  Eigen::MatrixXd visible_added = data.visible(Eigen::all, added_points);

  // Calculate the sum along columns and check if each sum is greater than zero
  Eigen::VectorXi affected_views = (visible_added.colwise().sum().array() > 0).cast<int>();
  std::cout<<"yessir"<<std::endl;
  throw std::exception();

}



std::pair<int, Eigen::VectorXi> FactorCompletion::try_adding_points(Eigen::VectorXi eligibles, int level_views){
  std::cout<<"try adding points entry"<<std::endl;
  Eigen::VectorXi added = Eigen::VectorXi::Zero(eligibles.size());
  int nums_added =0;
  auto& known_views = camera_variables.negative_pathway;

  std::cout<<"lets see 1"<<std::endl;
  for (int idx =0; idx<1; ++idx){
    std::vector<int> visible_views;
    for(int i=0; i<known_views.size(); i++){
      if(data.visible(i, eligibles[idx])){
        visible_views.push_back(i);
      }
    }
    Eigen::VectorXi visible_views_eig(visible_views.size());
    for(int i=0; i<visible_views.size(); i++){
      visible_views_eig(i) = visible_views[i];
    }
    Eigen::VectorXi inlier_views;
    Eigen::VectorXd best_estimate;
    if(Options::ROBUST_ESTIMATION){
      EstimatedRobustPoints* estim_views = new EstimatedRobustPoints(data, camera_variables,visible_views_eig, eligibles[idx], rejected_points(eligibles[idx]), level_views);
      best_estimate.resize(estim_views->best_estimate.size());
      inlier_views.resize(estim_views->best_inliers.size());
      best_estimate = estim_views->best_estimate;
      inlier_views= estim_views->best_inliers;
      std::cout<<"estim"<<std::endl;
      delete(estim_views);
      
      
    }else{
      EstimatedPoints* estim_views = new EstimatedPoints(data, camera_variables.cameras, visible_views_eig, eligibles[idx], visible_views_eig.size());
      best_estimate.resize(estim_views->estim.size());
      inlier_views.resize(visible_views.size());
      best_estimate = estim_views->estim;
      inlier_views= visible_views_eig;
      delete(estim_views);
    }

    if(best_estimate.size()==0){
      std::cout<<"if"<<std::endl;
      rejected_points(eligibles[idx]) = visible_views_eig.size();
    }
    else{
      std::cout<<"else"<<std::endl;
      camera_variables.points.col(eligibles[idx]) = best_estimate.block(0,0,4,1);
      inliers(inlier_views.array(),eligibles[idx]) = true;
      added(idx)=true;
      last_path++;
      nums_added++;
      camera_variables.pathway(last_path)= eligibles(idx);
      camera_variables.positive_pathway.conservativeResize(camera_variables.positive_pathway.size()+1);
      camera_variables.positive_pathway(camera_variables.positive_pathway.size()-1) = eligibles[idx];
      camera_variables.fixed.push_back(std::vector<int>(inlier_views.data(), inlier_views.data() + inlier_views.size()));
      rejected_points(eligibles[idx])= 0;
    }
  }
  return {nums_added, added};
}

int FactorCompletion::try_adding_views(std::vector<int> eligibles, int level_views){
  int nums_added = 0;
  auto& known_points = camera_variables.positive_pathway;
  for (int idx =0; idx<eligibles.size(); ++idx){

    std::vector<int> visible_pts;

    for(int i=0; i<known_points.size(); i++){
      if(data.visible(eligibles[idx],i)){
        visible_pts.push_back(i);
      }
    }
    Eigen::VectorXi visible_pts_eig(visible_pts.size());
    for(int i=0; i<visible_pts.size(); i++){
      visible_pts_eig(i) = visible_pts[i];
    }

    Eigen::VectorXi inlier_points;
    Eigen::VectorXd best_estimate;
    if(Options::ROBUST_ESTIMATION){
      EstimatedRobustViews estim_views = EstimatedRobustViews(data, camera_variables,visible_pts_eig, eligibles[idx], rejected_views(eligibles[idx]), level_views);
      best_estimate.resize(estim_views.best_estimate.size());
      inlier_points.resize(estim_views.best_inliers.size());
      best_estimate = estim_views.best_estimate;
      inlier_points = estim_views.best_inliers;
      // std::cout<<inlier_points<<std::endl;
      
      
    }else{
      EstimatedViews estim_views = EstimatedViews(data, camera_variables.points, visible_pts_eig, eligibles[idx], visible_pts_eig.size());
      best_estimate.resize(estim_views.estim.size());
      inlier_points.resize(visible_pts_eig.size());
      best_estimate = estim_views.estim;
      inlier_points = visible_pts_eig;
    }
    Eigen::Matrix<double, 3, 4> camera;
    if (best_estimate.size()==0) {
      rejected_views(eligibles[idx]) = visible_pts_eig.size();
    }
    else{
      camera << best_estimate.segment<4>(0).transpose(),
             best_estimate.segment<4>(4).transpose(),
             best_estimate.segment<4>(8).transpose();

      camera_variables.cameras.block<3, 4>(3 * eligibles[idx], 0) = camera;


      // std::cout<<camera_variables.cameras<<std::endl;
      inliers(eligibles[idx], inlier_points.array()) = true;
      nums_added++;
      last_path++;
      camera_variables.pathway(last_path) = -(eligibles[idx]+1);
      camera_variables.negative_pathway.conservativeResize(camera_variables.negative_pathway.size()+1);
      camera_variables.negative_pathway(camera_variables.negative_pathway.size()-1) = eligibles[idx];
      // std::cout<<"yup"<<std::endl;
      camera_variables.fixed.push_back(std::vector<int>(inlier_points.data(), inlier_points.data() + inlier_points.size()));
      rejected_views(eligibles[idx]) = 0;
    }
    
  }
  return nums_added;
}



Eigen::VectorXi FactorCompletion::search_eligible_points(int required_visible_eligibility, Eigen::VectorXi rejected_points, int iter) {

  int num_points = data.visible.cols();
  auto known_points = camera_variables.positive_pathway;
  auto known_views = camera_variables.negative_pathway;

  Eigen::VectorXi unknown_points = Eigen::VectorXi::Ones(num_points);
  for (int point : known_points) {
    unknown_points[point] = false;
  }
  
  Eigen::VectorXd num_visible_views = Eigen::VectorXd::Zero(num_points);
  for (int i = 0; i < num_points; ++i) {
    if (unknown_points[i]) {
      int count = 0;
      for (int view : known_views) {
        count = count + data.visible(view, i);
      }
      num_visible_views[i] = count;
    }
  }

  for(int i =0; i<num_visible_views.size(); i++){
    if(unknown_points(i)){
      num_visible_views[i]=num_visible_views[i];
    }
    else{
      num_visible_views[i]=false;
    }
  }

  Eigen::VectorXi rejected_unknown = rejected_points(unknown_points);
  Eigen::VectorXd::Scalar req_vis_elig_scalar = static_cast<Eigen::VectorXd::Scalar>(required_visible_eligibility);

  Eigen::VectorXd rejected_unknown_cast = rejected_unknown.cast<Eigen::VectorXd::Scalar>();
  Eigen::VectorXi condition1 = (num_visible_views.array() > rejected_unknown_cast.array()).cast<int>();
  for(int i = 0; i<num_visible_views.size(); i++){
    if (num_visible_views(i)>rejected_unknown(i)){
      condition1(i) = true;
    }
  }
  Eigen::VectorXi condition2 = (num_visible_views.array() >= req_vis_elig_scalar).cast<int>();
  Eigen::VectorXi eligible_unknown = condition1.cwiseProduct(condition2);

  // Eigen::VectorXi unknown(eligible_unknown.size());
  for(int i =0; i<unknown_points.size();i++){
    if(unknown_points(i)){
      unknown_points(i) = eligible_unknown[i];
    }
  }
  Eigen::VectorXi eligible_points((unknown_points.array()>0).count());
  int count =0;
  for(int i=0; i<unknown_points.size(); i++){
    if(unknown_points(i)){
      eligible_points(count)= i;
      count++;
    }
  }

  return eligible_points; 
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
  // deletion of heap allocated pvs scores
  // for(auto i: pvs_scores){
  //   delete(i);
  // }
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

  last_path = pathwayRef.size()-1;
  int new_zeros_length = num_views + num_points - last_path;

  pathwayRef.conservativeResize(Eigen::NoChange,
                                pathwayRef.cols() + new_zeros_length);
  pathwayRef.tail(new_zeros_length).setZero();

  for (int i = fixedRef.size() - new_zeros_length-1; i < fixedRef.size(); ++i) {
    fixedRef.push_back(std::vector<int>{});
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
