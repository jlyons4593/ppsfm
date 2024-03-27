#include "Helper.hpp"
#include "estimatedViews.h"
#include "dataStructures.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include "estimate_robust_views.h"
#include "options.h"

EstimatedRobustViews::EstimatedRobustViews(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables,Eigen::VectorXi known_points,int new_view, int num_rejected, int level):data(data), camera_variables(camera_variables){

    Eigen::Vector3i idx_view(3*new_view, 3*new_view+1, 3*new_view+2);
    double log_confidence = std::log(1-Options::CONFIDENCE/100);
    // std::cout<<log_confidence<<std::endl;
    double eps = std::numeric_limits<double>::epsilon();
    double best_score =  std::numeric_limits<double>::infinity();
    best_inliers = Eigen::VectorXi::Zero(known_points.size());
    int iteration_number = 0;
    int max_iterations = Options::MAX_ITERATION_ROBUST(1);


    //Main work loop
    // while (iteration_number<=0){
    while (iteration_number<=max_iterations){
        Eigen::VectorXi subset = Helper::random_subset(known_points, num_rejected, Options::MINIMAL_VIEW(level));

        iteration_number++;
        
        Eigen::VectorXi subset_points(subset.size());
        for(int i =0; i<subset.size(); i++){
            subset_points(i) = known_points(subset(i));
        }

        EstimatedViews* estim_views = new EstimatedViews(data, camera_variables.points, subset_points, new_view, 6);
        Eigen::VectorXd estim = estim_views->estim;
        Eigen::MatrixXd sys= estim_views->sys;
        delete(estim_views);
        

        Eigen::MatrixXd temp = sys*estim ;
        double threshold = temp.norm()/std::sqrt(sys.rows());
        if (estim.size()>0 && threshold<= Options::SYSTEM_THRESHOLD){
            Helper::InlierResults result = find_inliers(estim,idx_view, known_points);

            if (result.inliers.size() >= Options::MINIMAL_VIEW[level] && result.score< 7 * best_score) {
                Eigen::VectorXi known_points_subset(result.inliers.size());
                for(int i=0; i<result.inliers.size();i++){
                    known_points_subset(i) = known_points(result.inliers[i]);
                }

                EstimatedViews* estim_views2 = new EstimatedViews(data, camera_variables.points, known_points_subset, new_view, result.inliers.size());

                Eigen::VectorXd estim2 = estim_views2->estim;
                Eigen::MatrixXd sys2 = estim_views2->sys;
                delete(estim_views2);

                Eigen::MatrixXd temp2 = sys2*estim2 ;
                double threshold2 = temp2.norm()/std::sqrt(sys2.rows());

                if (estim2.size()>0 && threshold2<= Options::SYSTEM_THRESHOLD){
                    double score = computeScore(estim2, idx_view, known_points, result.inliers);
                    if(score<best_score){
                        best_score = score;
                        best_estimate = estim2;
                        best_inliers = known_points(result.inliers);


                        double ratio = static_cast<double>(result.inliers.size()) / known_points.size();
                        double prob = std::max(eps, std::min(1 - eps, 1 - std::pow(ratio,Options::MINIMAL_VIEW[level])));

                        max_iterations= std::min(static_cast<int>(std::ceil((log_confidence) / std::log(prob))), Options::MAX_ITERATION_ROBUST[1]);
                    }
                }
            }
        }
        iteration_number= iteration_number+1;
    }
    // std::cout<<"best_inliers:"<<std::endl;
    // std::cout<<best_inliers<<std::endl;
    // std::cout<<"best_estimate:"<<std::endl;
    // std::cout<<best_estimate<<std::endl;
    std::cout<<"best_score:"<<std::endl;
    std::cout<<best_score<<std::endl;
}
double EstimatedRobustViews::computeScore(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points, std::vector<int> inliers){
    std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, idx_view, known_points);

    double score = reproj.first(inliers).sum() + (known_points.size() - inliers.size()) * Options::OUTLIER_THRESHOLD;
    // std::cout<<score<<std::endl;
    return score;
    
}

std::pair<Eigen::VectorXd, Eigen::RowVectorXd> EstimatedRobustViews::compute_reproj(Eigen::VectorXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points)
{    

    // Denormalize the estimated camera
    Eigen::MatrixXd A(idx_view.size(), data.normalisations.cols());
    for(int i = 0; i < idx_view.size(); ++i) {
        A.row(i) = data.normalisations.row(idx_view(i));
    }
    Eigen::MatrixXd temp(3, 4);
    temp << estimation.segment(0, 4).transpose(),
         estimation.segment(4, 4).transpose(),
         estimation.segment(8, 4).transpose();

    Eigen::MatrixXd camera = A.inverse() * temp;
    // std::cout<<"test"<<std::endl;
    Eigen::MatrixXd scaledMeasurements = camera * camera_variables.points(Eigen::all, known_points);
    Eigen::RowVectorXd depths = scaledMeasurements.row(2);

    Eigen::MatrixXd reprojections = scaledMeasurements.cwiseQuotient(depths.replicate(scaledMeasurements.rows(), 1));

    Eigen::MatrixXd img_meas_subset = data.image_measurements(idx_view, Eigen::all);
    Eigen::MatrixXd img_meas_subset_points = img_meas_subset(Eigen::all, known_points);

    Eigen::MatrixXd reproj_err = reprojections - img_meas_subset_points;
    Eigen::VectorXd reproj_err_squared = reproj_err.row(0).array().square() + reproj_err.row(1).array().square();
    Eigen::VectorXd reproj_err_final = reproj_err_squared.array().sqrt();

    return {reproj_err_final,depths};
}

Helper::InlierResults EstimatedRobustViews::find_inliers(Eigen::MatrixXd estimation,Eigen::Vector3i idx_view, Eigen::VectorXi known_points){

    Helper::InlierResults results;

    std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, idx_view, known_points);

    Eigen::VectorXd reproj_errs = reproj.first;
    Eigen::VectorXd inliers = (reproj_errs.array() <= Options::OUTLIER_THRESHOLD).cast<double>();
    std::vector<int> inliers_indices;

    for (int i = 0; i < inliers.size(); ++i) {
        if (inliers[i] == 1) {
            inliers_indices.push_back(i);
        }
    }

    results.score = reproj_errs(inliers_indices).sum() + (known_points.size()-  inliers_indices.size()) * Options::OUTLIER_THRESHOLD;

    results.inliers = inliers_indices;
    return results;
}

