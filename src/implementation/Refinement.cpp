#include "Refinement.h"
#include "estimatedViews.h"
#include "estimatedPoints.h"
#include "dataStructures.hpp"
#include "options.h"
#include <iostream>
#include <limits>
#include <vector>

Refinement::Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> visible, int init_refine, int last_path, bool start_cameras, int type)
    :data(data), camera_variables(camera), visible(visible){

        
    Eigen::VectorXi pathway_segment = camera_variables.pathway.segment(init_refine-1, last_path+1 );

    std::vector<std::vector<int>> fixed;
    fixed.assign(std::next(camera_variables.fixed.begin(), init_refine - 1), std::next(camera_variables.fixed.begin(), last_path+1));

    int max_iter;
    double min_change;
    if(type ==2){
        max_iter = Options::MAX_ITER_FINAL_REFINE;   
        min_change = Options::MIN_CHANGE_FINAL_REFINE;
    }
    else if (type==1){
        max_iter = Options::MAX_ITERATION_REFINEMENT;
        min_change = Options::MIN_CHANGE_GLOBAL_REFINE;
    }
    else{
        max_iter = Options::MAX_ITERATION_REFINEMENT;
        min_change = Options::MIN_CHANGE_LOCAL_REFINEMENT;
    }
    Eigen::VectorXi idx_points(camera_variables.positive_pathway.size());  
    // std::cout<<camera_variables.negative_pathway<<std::endl;
    Eigen::VectorXi idx_cameras(camera_variables.negative_pathway.size()); 
    int camera_count =0;
    int point_count =0;
    for(int i =0; i<pathway_segment.size();i++){
        if(pathway_segment(i)<0){
            idx_cameras(camera_count) = i;
            camera_count++;
        }else{
            idx_points(point_count) = i;
            point_count++;
        }
    }

    Eigen::VectorXi selected_pathway = pathway_segment(idx_cameras);
    Eigen::VectorXi negated_scaled_pathway = selected_pathway * -3;

    Eigen::MatrixXi replicated_pathway = negated_scaled_pathway.replicate(1, 3);
    replicated_pathway = replicated_pathway.array()-1;
    replicated_pathway.resizeLike(Eigen::MatrixXi::Ones(idx_cameras.size(), 3));

    Eigen::RowVectorXi sequence(3);
    sequence << 2, 1, 0;

    Eigen::MatrixXi replicated_sequence = sequence.replicate(idx_cameras.size(), 1);
    Eigen::VectorXi known_cams3 = (replicated_pathway - replicated_sequence).reshaped();

    std::vector<int> rowSize;
    for (int i =0; i<visible.cols(); i++){
        rowSize.push_back(i);
    }
    Eigen::VectorXi positive_path = camera_variables.positive_pathway;

    std::sort(rowSize.begin(), rowSize.end());
    std::sort(positive_path.begin(), positive_path.end());

    // Compute difference for visible rows to be changed
    std::vector<int> diff;
    std::set_difference(rowSize.begin(), rowSize.end(), positive_path.begin(), positive_path.end(), std::back_inserter(diff));

    for(int i=0; i<diff.size(); i++){
        visible.col(diff[i]).array() = 0;
    }

    // Compute difference for visible cols to be changed
    std::vector<int> colSize;
    for (int i =0; i<visible.rows(); i++){
        colSize.push_back(i); 
    }


    std::vector<int> diff2;
    std::set_difference(colSize.begin(), colSize.end(), camera_variables.negative_pathway.begin(), camera_variables.negative_pathway.end(), std::back_inserter(diff2));
    for(int i=0; i<diff2.size(); i++){
        visible.row(diff2[i]).array() = 0;
    }

    // While loop vars

    int iteration_number=0;
    double change_cameras =std::numeric_limits<double>::infinity();
    double change_points =std::numeric_limits<double>::infinity();
    // std::cout<<positive_path<<std::endl;

    // Optional timer later

    while(iteration_number<max_iter){
    // while(iteration_number<1){
        Eigen::MatrixXd old_cameras(known_cams3.size(), camera_variables.cameras.cols()); 
        for(int i=0; i<known_cams3.size(); i++){
            old_cameras.row(i) = camera_variables.cameras.row(known_cams3(i));
        }
        
        Eigen::MatrixXd old_points(camera_variables.points.rows(), positive_path.size()); 
        for(int i=0; i<camera_variables.points.rows(); i++){
            old_points.row(i) = camera_variables.points.block(i,0,1,positive_path.size());
        }

        if(start_cameras){
            reestimate_all_views(pathway_segment, idx_cameras);
            reestimate_all_points(pathway_segment, idx_points);
        }
        else{
            reestimate_all_points(pathway_segment, idx_points);
            reestimate_all_views(pathway_segment, idx_cameras);
        }
        
        // Checking cams
        Eigen::MatrixXd diff = old_cameras - camera_variables.cameras(known_cams3, Eigen::all);
        double norm_diff = diff.norm();
        int numel_old_cameras = old_cameras.size();
        double sqrt_numel = std::sqrt(numel_old_cameras);
        double change_cam = norm_diff / sqrt_numel;
        // Checking points
        Eigen::MatrixXd diff_points = old_points - camera_variables.points(Eigen::all, positive_path);
        double norm_diff_points = diff_points.norm();
        int numel_old_points = old_points.size();
        double sqrt_numel_points = std::sqrt(numel_old_points);
        double change_pts = norm_diff_points / sqrt_numel_points;

        if(iteration_number>Options::MAX_ITERATION_REFINEMENT&& (change_pts+change_cam)/2<min_change){
            break;

        }
        // std::cout<<change_pts<<std::endl;
        iteration_number++;
    }

}
    
void Refinement::reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points){

    for (int i =0; i<idx_points.size(); i++){
    // for (int i =0; i<1; i++){
        int idx = pathway(idx_points(i))-1;
        std::vector<int> fixed_views;
        fixed_views = camera_variables.fixed[idx_points[i]];

        Eigen::VectorXi visible_views = visible.col(idx).cast<int>();
        // std::cout<<visible_views<<std::endl;
        for(auto i: fixed_views){
            visible_views(i) = false;
        }
        // std::cout<<visible_views<<std::endl;


        Eigen::VectorXi fixed_views_eigen = Eigen::Map<Eigen::VectorXi>(fixed_views.data(), fixed_views.size());
        
        std::vector<int> test;
        for(int i=0 ; i <visible_views.size(); i++){
            if(visible_views(i)){
                test.push_back(i);
            }
        }
        Eigen::VectorXi visible_views2 = Eigen::Map<Eigen::VectorXi>(test.data(), test.size());

        fixed_views_eigen.conservativeResize(fixed_views_eigen.size() + visible_views2.size());
        fixed_views_eigen.tail(visible_views2.size()) = visible_views2;

        

        EstimatedPoints* estim = new EstimatedPoints(data,camera_variables.cameras, fixed_views_eigen, idx, fixed_views.size());
        Eigen::VectorXd estimation = estim->estim;
        // std::cout<<estimation<<std::endl;
        delete(estim);

        if(estimation.size()==0){
            std::cout<<"Estimation point is empty"<<std::endl;
        }
        else{
            assert(estimation.size()==4);
            camera_variables.points.col(idx) = estimation;
        }
    }
}

void Refinement::reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras){

    // for (int i =0; i<idx_cameras.size(); i++){
    for (int i =0; i<1; i++){
        int idx = (-pathway(idx_cameras(i)))-1;
        std::vector<int> fixed_points;
        fixed_points= camera_variables.fixed[idx_cameras[i]];
        Eigen::VectorXi visible_points= visible.row(idx).cast<int>();
        for(auto i: fixed_points){
            visible_points(i) = false;
        }

        Eigen::VectorXi fixed_points_eigen = Eigen::Map<Eigen::VectorXi>(fixed_points.data(), fixed_points.size());
        
        std::vector<int> test;
        for(int i=0 ; i <visible_points.size(); i++){
            if(visible_points(i)){
                test.push_back(i);
            }
        }

        Eigen::VectorXi visible_points2 = Eigen::Map<Eigen::VectorXi>(test.data(), test.size());
        fixed_points_eigen.conservativeResize(fixed_points_eigen.size() + visible_points2.size());
        fixed_points_eigen.tail(visible_points2.size()) = visible_points2;

        EstimatedViews* estim = new EstimatedViews(data, camera_variables.points, fixed_points_eigen, idx , fixed_points.size());
        // std::cout<<estim->estim<<std::endl;
        Eigen::VectorXd estimation= estim->estim;
        delete(estim);

        if(estimation.size()==0){
            std::cout<<"Estimation view is empty"<<std::endl;
        }
        else{
            // std::cout<<estimation<<std::endl;
            assert(estimation.size()==12);
            Eigen::MatrixXd mat = estimation;
            mat.resize(4,3);

            camera_variables.cameras.block(idx*3,0,3,4) = mat.transpose();

        }
        
    }

}

