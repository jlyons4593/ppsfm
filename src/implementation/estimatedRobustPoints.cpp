#include "estimatedRobustPoints.h"
#include "dataStructures.hpp"
#include "Helper.hpp"
#include "estimatedPoints.h"
#include "estimatedViews.h"
#include "dataStructures.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include "options.h"

EstimatedRobustPoints::EstimatedRobustPoints(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables, Eigen::VectorXi& known_views, int new_point, int num_rejected, int level): data(data), camera_variables(camera_variables){

    
    std::sort(known_views.data(), known_views.data() + known_views.size());
    int num_known_views = known_views.size();
    Eigen::MatrixXd denorm_cams(3 * num_known_views, 4);
    denorm_cams.setConstant(std::numeric_limits<double>::quiet_NaN());
    // if (new_point == 258){
    //     std::cout<<"yup0"<<std::endl;
    // }

    // have a look at denorm cam normalisation block can be the same

    // Eigen::MatrixXd cameras_block = camera_variables.cameras.block(19*3,0,3,camera_variables.cameras.cols());
        
    for(int j=0; j<known_views.size(); j++){
        Eigen::MatrixXd normalisation_block = data.normalisations.block(known_views(j)*3,0,3,data.normalisations.cols());
        Eigen::MatrixXd cameras_block = camera_variables.cameras.block(known_views(j)*3,0,3,camera_variables.cameras.cols());
        denorm_cams.block(3 * j, 0, 3, 4) = normalisation_block.colPivHouseholderQr().solve(cameras_block);
    }
    // if (new_point == 258){
    //     std::cout<<"yup1"<<std::endl;
    // }


    Eigen::VectorXi idx_views(3 * num_known_views);
    for (int i = 0; i < num_known_views; ++i) {
        idx_views.segment(3 * i, 3).setConstant(3 * known_views[i]);
    }

    // if (new_point == 258){
    //     std::cout<<"yup2"<<std::endl;
    // }

    for (int i = 0; i < idx_views.size(); i += 3) {
        idx_views[i] += 1; 
        if (i + 1 < idx_views.size()) {
            idx_views[i + 2] += 2; 
        }
    }

    // if (new_point == 258){
    //     std::cout<<"yup3"<<std::endl;
    // }
    double log_confidence = std::log(1-Options::CONFIDENCE/100);
    // std::cout<<log_confidence<<std::endl;
    double eps = std::numeric_limits<double>::epsilon();
    double best_score =  std::numeric_limits<double>::infinity();
    best_inliers = Eigen::VectorXi::Zero(known_views.size());
    int iteration_number = 0;
    // int max_iterations = Options::MAX_ITERATION_ROBUST(1);

    int max_iterations = 20;
    // if (new_point == 258){
    //     std::cout<<"yup4"<<std::endl;
    // }

    //Main work loop
    // std::cout<<"iteration number"<<std::endl;
    while (iteration_number<=max_iterations){
        // std::cout<< iteration_number<< " ";
        Eigen::VectorXi subset = Helper::random_subset(known_views, num_rejected, Options::MINIMAL_POINT(level));
        // if (new_point == 258){
        //     std::cout<<"subset"<<std::endl;
        // std::cout<<subset<<std::endl;
        // }
        
        Eigen::VectorXi subset_views(subset.size());
        for(int i =0; i<subset.size(); i++){
            subset_views(i) = known_views(subset(i));
        }
        // if (new_point == 258){
        //     std::cout<<"yup4"<<std::endl;
        //     std::cout<<subset_views<<std::endl;
        // }

        EstimatedPoints estim_points= EstimatedPoints(data, camera_variables.cameras, subset_views, new_point, 2);
        Eigen::VectorXd estim = estim_points.estim;
    // if (new_point == 258){
    //     std::cout<<"yup5"<<std::endl;
    //     std::cout<<estim<<std::endl;
    // }

        Eigen::MatrixXd sys= estim_points.sys;
        Eigen::MatrixXd temp = sys*estim ;
        double threshold = temp.norm()/std::sqrt(sys.rows());


        // std::cout<<"pre if"<<std::endl;
        if (estim.size()>0 && threshold<= Options::SYSTEM_THRESHOLD){
            Helper::InlierResults result = find_inliers(estim ,denorm_cams, idx_views, new_point);


            if (result.inliers.size() >= Options::MINIMAL_POINT[level] && result.score< 7 * best_score) {
                // if (new_point == 258){
                //     std::cout<<"yup6"<<std::endl;
                // }

                // std::cout<<"yup1"<<std::endl;
                Eigen::VectorXi known_views_subset(result.inliers.size());
                for(int i=0; i<result.inliers.size();i++){
                    known_views_subset(i) = known_views(result.inliers[i]);
                }

                EstimatedPoints estim_views2 = EstimatedPoints(data, camera_variables.cameras, known_views_subset, new_point, result.inliers.size());
                Eigen::VectorXd estim2 = estim_views2.estim;
                Eigen::MatrixXd sys2 = estim_views2.sys;

                Eigen::MatrixXd temp2 = sys2*estim2 ;
                double threshold2 = temp2.norm()/std::sqrt(sys2.rows());

                if (estim2.size()>0 && threshold2<= Options::SYSTEM_THRESHOLD){
                    double score = computeScore(estim2, denorm_cams,  idx_views, new_point, result.inliers);
                    if(score<best_score){
                        best_score = score;
                        best_estimate = estim2;
                        best_inliers = known_views(result.inliers);

                        double ratio = static_cast<double>(result.inliers.size()) / known_views.size();
                        double prob = std::max(eps, std::min(1 - eps, 1 - std::pow(ratio,Options::MINIMAL_POINT[level])));

                        // max_iterations= std::min(static_cast<int>(std::ceil((log_confidence) / std::log(prob))), Options::MAX_ITERATION_ROBUST[0]);
                        max_iterations= std::min(static_cast<int>(std::ceil((log_confidence) / std::log(prob))),20);
                        // if(static_cast<int>(std::ceil((log_confidence) / std::log(prob)))){
                        //     std::cout<<"max iterations: "<<max_iterations<<std::endl;
                        // }
                    }
                }
            }
        }
        
        iteration_number++;
    }
    // std::cout<<std::endl;
    // std::cout<<"iteration number: "<<iteration_number<<std::endl;
        // if (new_point == 258){
        //     std::cout<<"best estim"<<std::endl;
        //     std::cout<<best_estimate<<std::endl;
        // }

}

double EstimatedRobustPoints::computeScore(Eigen::VectorXd& estimation, Eigen::MatrixXd& cams, Eigen::VectorXi& idx_view, int new_point, std::vector<int> inliers){
    std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, idx_view,cams, new_point);
    Eigen::VectorXd reproj_errs = reproj.first;

    double score = reproj_errs(inliers).sum() + ((idx_view.size()/3.00) - inliers.size())*Options::OUTLIER_THRESHOLD;
    return score;
    
}

std::pair<Eigen::VectorXd, Eigen::RowVectorXd> EstimatedRobustPoints::compute_reproj(Eigen::VectorXd& estimation,Eigen::VectorXi& idx_view, Eigen::MatrixXd& cameras, int new_point)
{    
    // could be issues with scaled measurements/ cameras and estimations being incorrect

    std::sort(idx_view.data(), idx_view.data() + idx_view.size());
    Eigen::MatrixXd scaled_measurements = cameras * estimation;
    int depth_size= scaled_measurements.size()  / 3;

// Create a new vector to store every 3rd element
    Eigen::VectorXd depths(depth_size);

    // Copy every 3rd element from the original vector to the new vector
    for (int i = 1; i <= depth_size; ++i) {
        depths(i-1) = scaled_measurements((3 * i)-1);
    }
    // std::cout<<"depths"<<std::endl;
    // std::cout<<depths<<std::endl;
    Eigen::VectorXd kronResult(scaled_measurements.size()) ;
    for(int i=0; i<depths.size(); i++){
        kronResult(i*3) = depths(i);
        kronResult(i*3+1) = depths(i);
        kronResult(i*3+2) = depths(i);
    }
    
    Eigen::VectorXd reprojections = scaled_measurements.array() / kronResult.array();
    Eigen::VectorXd img_meas_subset_points(idx_view.size());
    for(int i = 0 ; i<idx_view.size();i++){
        img_meas_subset_points(i) = data.image_measurements(idx_view(i), new_point);
    }
    Eigen::MatrixXd reproj_err = reprojections - img_meas_subset_points;
    // std::cout<<"reproj err"<<std::endl;
    // std::cout<<reproj_err<<std::endl;
    int newSize = reproj_err.size() / 3;
    Eigen::VectorXd result(newSize);


    for (int i = 0; i < newSize; ++i) {
        double val1 = reproj_err(3 * i);     // Picks 1st, 4th, 7th, etc.
        double val2 = reproj_err(3 * i + 1); // Picks 2nd, 5th, 8th, etc.
        result(i) = std::sqrt(val1 * val1 + val2 * val2);
    }
    return {result,depths};
}

Helper::InlierResults EstimatedRobustPoints::find_inliers(Eigen::VectorXd& estimation,Eigen::MatrixXd& cameras,Eigen::VectorXi idx_view, int new_point){

    Helper::InlierResults results;
    std::pair<Eigen::VectorXd, Eigen::RowVectorXd> reproj = compute_reproj(estimation, idx_view, cameras, new_point);
    Eigen::VectorXd reproj_errs = reproj.first;
    // std::cout<<"reprojection errors"<<std::endl;
    // std::cout<<reproj_errs<<std::endl;
    Eigen::VectorXd inliers = (reproj_errs.array() <= Options::OUTLIER_THRESHOLD).cast<double>();
    // std::cout<<inliers<<std::endl;
    std::vector<int> inliers_indices;

    for (int i = 0; i < inliers.size(); ++i) {
        if (inliers[i] == 1) {
            inliers_indices.push_back(i);
        }
    }

    results.score = reproj_errs(inliers).sum() + ((idx_view.size()/3.00) - inliers.size())*Options::OUTLIER_THRESHOLD;
    results.inliers = inliers_indices;
    results.reproj_errors= reproj_errs(inliers);
    return results;
}

