#include "estimatedRobustPoints.h"
#include "dataStructures.hpp"
#include "Helper.hpp"
#include "estimatedPoints.h"
#include "estimatedViews.h"
#include "dataStructures.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include "options.h"

EstimatedRobustPoints::EstimatedRobustPoints(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables, Eigen::VectorXi& known_views, int new_point, int num_rejected, int level){

    
    int num_known_views = known_views.size();
    Eigen::MatrixXd denorm_cams(3 * num_known_views, 4);
    denorm_cams.setConstant(std::numeric_limits<double>::quiet_NaN());

    //Translate this
    //

    Eigen::VectorXi idx_views(3 * num_known_views);
    for (int i = 0; i < num_known_views; ++i) {
        idx_views.segment(3 * i, 3).setConstant(3 * known_views[i]);
    }


    for (int i = 0; i < idx_views.size(); i += 3) {
        idx_views[i] += 1; 
        if (i + 1 < idx_views.size()) {
            idx_views[i + 2] += 2; 
        }
    }

    double log_confidence = std::log(1-Options::CONFIDENCE/100);
    // std::cout<<log_confidence<<std::endl;
    double eps = std::numeric_limits<double>::epsilon();
    double best_score =  std::numeric_limits<double>::infinity();
    best_inliers = Eigen::VectorXi::Zero(known_views.size());
    int iteration_number = 0;
    int max_iterations = Options::MAX_ITERATION_ROBUST(1);


    //Main work loop
    // while (iteration_number<=0){
    while (iteration_number<=max_iterations){
        Eigen::VectorXi subset = Helper::random_subset(known_views, num_rejected, Options::MINIMAL_VIEW(level));

        
        Eigen::VectorXi subset_views(subset.size());
        for(int i =0; i<subset.size(); i++){
            subset_views(i) = known_views(subset(i));
        }

        EstimatedPoints* estim_points= new EstimatedPoints(data, camera_variables.cameras, subset_views, new_point, 2);
        Eigen::VectorXd estim = estim_points->estim;
        Eigen::MatrixXd sys= estim_points->sys;
        delete(estim_points);
        Eigen::MatrixXd temp = sys*estim ;
        double threshold = temp.norm()/std::sqrt(sys.rows());

        if (estim.size()>0 && threshold<= Options::SYSTEM_THRESHOLD){
            // Helper::InlierResults result = find_inliers(estim, idx_views, denorm_cams);
            //
            // if (result.inliers.size() >= Options::MINIMAL_VIEW[level] && result.score< 7 * best_score) {
            //     Eigen::VectorXi known_views_subset(result.inliers.size());
            //     for(int i=0; i<result.inliers.size();i++){
            //         known_views_subset(i) = known_views(result.inliers[i]);
            //     }
            //
            //     EstimatedViews* estim_views2 = new EstimatedViews(data, camera_variables.points, known_views_subset, new_point, result.inliers.size());
            //     Eigen::VectorXd estim2 = estim_views2->estim;
            //     Eigen::MatrixXd sys2 = estim_views2->sys;
            //     delete(estim_views2);
            //
            //     Eigen::MatrixXd temp2 = sys2*estim2 ;
            //     double threshold2 = temp2.norm()/std::sqrt(sys2.rows());
            //
            //     if (estim2.size()>0 && threshold2<= Options::SYSTEM_THRESHOLD){
            //         double score = computeScore(estim2, denorm_cams,  idx_views, known_views, result.inliers);
            //         if(score<best_score){
            //             best_score = score;
            //             best_estimate = estim2;
            //             best_inliers = known_views(result.inliers);
            //
            //             double ratio = static_cast<double>(result.inliers.size()) / known_views.size();
            //             double prob = std::max(eps, std::min(1 - eps, 1 - std::pow(ratio,Options::MINIMAL_VIEW[level])));
            //
            //             max_iterations= std::min(static_cast<int>(std::ceil((log_confidence) / std::log(prob))), Options::MAX_ITERATION_ROBUST[1]);
            //         }
            //     }
            // }
        }
        iteration_number++;
    }
}

// double EstimatedRobustPoints::computeScore(Eigen::VectorXd& estimation, Eigen::MatrixXd& cams, Eigen::VectorXi& idx_view, Eigen::VectorXi& known_views, std::vector<int> inliers){
    // std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, cams, idx_view, known_views);
    // double score = reproj.first(inliers).sum() + (known_views.size() - inliers.size()) * Options::OUTLIER_THRESHOLD;

    // return score;
    
// }

// std::pair<Eigen::VectorXd, Eigen::RowVectorXd> EstimatedRobustViews::compute_reproj(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points)
// {    
//
//     // Denormalize the estimated camera
//     Eigen::MatrixXd A(idx_view.size(), data.normalisations.cols());
//     for(int i = 0; i < idx_view.size(); ++i) {
//         A.row(i) = data.normalisations.row(idx_view(i));
//     }
//     Eigen::MatrixXd temp(3, 4);
//     temp << estimation.segment(0, 4).transpose(),
//          estimation.segment(4, 4).transpose(),
//          estimation.segment(8, 4).transpose();
//
//     Eigen::MatrixXd camera = A.inverse() * temp;
//     // std::cout<<"test"<<std::endl;
//     Eigen::MatrixXd scaledMeasurements = camera * camera_variables.points(Eigen::all, known_points);
//     Eigen::RowVectorXd depths = scaledMeasurements.row(2);
//
//     Eigen::MatrixXd reprojections = scaledMeasurements.cwiseQuotient(depths.replicate(scaledMeasurements.rows(), 1));
//
//     Eigen::MatrixXd img_meas_subset = data.image_measurements(idx_view, Eigen::all);
//     Eigen::MatrixXd img_meas_subset_points = img_meas_subset(Eigen::all, known_points);
//
//     Eigen::MatrixXd reproj_err = reprojections - img_meas_subset_points;
//     Eigen::VectorXd reproj_err_squared = reproj_err.row(0).array().square() + reproj_err.row(1).array().square();
//     Eigen::VectorXd reproj_err_final = reproj_err_squared.array().sqrt();
//
//     return {reproj_err_final,depths};
// }
//
// Helper::InlierResults EstimatedRobustViews::find_inliers(Eigen::MatrixXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points){
//
//     Helper::InlierResults results;
//
//     std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, idx_view, known_points);
//
//     Eigen::VectorXd reproj_errs = reproj.first;
//     Eigen::VectorXd inliers = (reproj_errs.array() <= Options::OUTLIER_THRESHOLD).cast<double>();
//     std::vector<int> inliers_indices;
//
//     for (int i = 0; i < inliers.size(); ++i) {
//         if (inliers[i] == 1) {
//             inliers_indices.push_back(i);
//         }
//     }
//
//     results.score = reproj_errs(inliers_indices).sum() + (known_points.size()-  inliers_indices.size()) * Options::OUTLIER_THRESHOLD;
//
//     results.inliers = inliers_indices;
//     return results;
// }

