#include "Refinement.h"
#include <chrono>
#include <exception>
#include <stdexcept>
#include "estimatedViews.h"
#include "estimatedPoints.h"
#include "dataStructures.hpp"
#include "options.h"
#include <iostream>
#include <limits>
#include <vector>

Eigen::MatrixXd Refinement::getCameras(){
    return camera_variables.cameras;
}
Eigen::MatrixXd Refinement::getPoints(){
    return camera_variables.points;
}
Refinement::Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> visible,Eigen::VectorXi pathway_segment,std::vector<std::vector<int>> new_fixed,    bool start_cameras, int type)
    :data(data), camera_variables(camera), visible(visible),pathway_segment(pathway_segment){

        std::cout<<"yup0"<<std::endl;
        fixed = new_fixed;
        // camera_variables.pathway.segment(init_refine,last_path - init_refine)

        std::cout<<"yup1"<<std::endl;
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
    std::cout<<"yup2"<<std::endl;
    Eigen::VectorXi idx_points((pathway_segment.array()>0).count());  
    // std::cout<<"test"<<std::endl;
    // std::cout<<"pathway segment"<<std::endl;
    // std::cout<<pathway_segment.array()<<std::endl;
    Eigen::VectorXi idx_cameras((pathway_segment.array()<0).count()); 
        std::cout<<"yup3"<<std::endl;
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
        // std::cout<<i<<std::endl;
    }
        std::cout<<"yup4"<<std::endl;
    // std::cout<<"yessir1"<<std::endl;

    Eigen::VectorXi selected_pathway = pathway_segment(idx_cameras);
    Eigen::VectorXi negated_scaled_pathway = selected_pathway * -3;

    Eigen::MatrixXi replicated_pathway = negated_scaled_pathway.replicate(1, 3);
    replicated_pathway = replicated_pathway.array()-1;
    replicated_pathway.resizeLike(Eigen::MatrixXi::Ones(idx_cameras.size(), 3));

    Eigen::RowVectorXi sequence(3);
    sequence << 2, 1, 0;

    Eigen::MatrixXi replicated_sequence = sequence.replicate(idx_cameras.size(), 1);
    Eigen::VectorXi known_cams3 = (replicated_pathway - replicated_sequence).reshaped();
        std::cout<<"yup5"<<std::endl;

    std::vector<int> rowSize;
    for (int i =0; i<visible.cols(); i++){
        rowSize.push_back(i);
    }

    std::vector<int> positiveValues;
    std::copy_if(pathway_segment.data(), pathway_segment.data() + pathway_segment.size(), std::back_inserter(positiveValues),
                 [](int value) { return value > 0; });
    // Now create a new Eigen vector from the std::vector of positive values
    Eigen::Map<Eigen::VectorXi> positive_path(positiveValues.data(), positiveValues.size());
    positive_path =positive_path.array()-1;

    std::vector<int> negativeValues;
    std::copy_if(pathway_segment.data(), pathway_segment.data() + pathway_segment.size(), std::back_inserter(negativeValues),
                 [](int value) { return value < 0; });

        std::cout<<"yup6"<<std::endl;
    // Now create a new Eigen vector from the std::vector of positive values
    Eigen::Map<Eigen::VectorXi> negative_path(negativeValues.data(), negativeValues.size());
    negative_path= -(negative_path.array()) -1;

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
    std::set_difference(colSize.begin(), colSize.end(), negative_path.begin(), negative_path.end(), std::back_inserter(diff2));
    for(int i=0; i<diff2.size(); i++){
        visible.row(diff2[i]).array() = 0;
    }

    // While loop vars

    int iteration_number=1;
    double change_cameras =std::numeric_limits<double>::infinity();
    double change_points =std::numeric_limits<double>::infinity();

    // Optional timer later
    auto start = std::chrono::high_resolution_clock::now();

        // std::cout<<"yup7"<<std::endl;

    // DEBUG WITH ITERATION NUMBER 1 the change camera vars and points vars are far too far away
    while(iteration_number<max_iter){
        std::cout<<"refinement iteration: "<<iteration_number<<std::endl;
    // while(iteration_number<1){
        Eigen::MatrixXd old_cameras(known_cams3.size(), camera_variables.cameras.cols()); 
        for(int i=0; i<known_cams3.size(); i++){
            old_cameras.row(i) = camera_variables.cameras.row(known_cams3(i));
        }
        
        // std::cout<<"yup8"<<std::endl;
        // Eigen::MatrixXd old_points(camera_variables.points.rows(), positive_path.size()); 
        // for(int i=0; i<camera_variables.points.rows(); i++){
        //     old_points.row(i) = camera_variables.points.block(i,0,1,positive_path.size());
        // }
        // std::cout<<positive_path.transpose()<<std::endl;
        Eigen::MatrixXd old_points = camera_variables.points(Eigen::all, positive_path);

        // std::cout<<"yup9"<<std::endl;

        if(start_cameras){
        // std::cout<<"reestimate views"<<std::endl;
            reestimate_all_views(pathway_segment, idx_cameras, fixed);
        // std::cout<<"reestimate points"<<std::endl;
            reestimate_all_points(pathway_segment, idx_points, fixed);
        }
        else{
            // std::cout<<"idx points"<<std::endl;
            // std::cout<<idx_points.transpose()<<std::endl;
            // std::cout<<"idx_cameras"<<std::endl;
            // std::cout<<idx_cameras.transpose()<<std::endl;
        // std::cout<<"reestimate points"<<std::endl;
            reestimate_all_points(pathway_segment, idx_points, fixed);
        // std::cout<<"reestimate views"<<std::endl;
            reestimate_all_views(pathway_segment, idx_cameras, fixed);
        }
        
        // Checking cams
        Eigen::MatrixXd diff = old_cameras - camera_variables.cameras(known_cams3, Eigen::all);
        double norm_diff = diff.norm();
        int numel_old_cameras = old_cameras.size();
        double sqrt_numel = std::sqrt(numel_old_cameras);
        change_cameras= norm_diff / sqrt_numel;
        // Checking points
        Eigen::MatrixXd diff_points = old_points - camera_variables.points(Eigen::all, positive_path);
        double norm_diff_points = diff_points.norm();
        int numel_old_points = old_points.size();
        double sqrt_numel_points = std::sqrt(numel_old_points);
        change_points= norm_diff_points / sqrt_numel_points;

        if(iteration_number>Options::MIN_ITERATION_REFINEMENT-1 && (change_points+change_cameras)/2<min_change){
            std::cout<<"broke out on iter "<<iteration_number<<std::endl;
            break;

        }
        // if(iteration_number>2){
        //     std::cout<<"greater than 3"<<std::endl;
        //     // throw std::exception();
        // }
        iteration_number++;
    }
    // std::cout<<"refinement "<<iteration_number <<" iterations' "<<std::endl;

    // if(Options::LOG_TIMER){
    //     auto end = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "Refinement execution time: " << duration.count() << " Milliseconds" << std::endl;
        // if(duration.count()>1000){
        //     std::cout<<"Long process"<<std::endl;
        // }
    // }
}
    
void Refinement::reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points, std::vector<std::vector<int>> fixed){

    for (int i =0; i<idx_points.size(); i++){
        auto start= std::chrono::high_resolution_clock::now();
        int idx = pathway(idx_points(i))-1;
        std::vector<int> fixed_views;
        fixed_views = fixed[idx_points[i]];
        Eigen::VectorXi visible_views = visible.col(idx).cast<int>();
        for(auto i: fixed_views){
            visible_views(i) = false;
        }

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



        EstimatedPoints estim = EstimatedPoints(data,camera_variables.cameras, fixed_views_eigen, idx, fixed_views.size());
        Eigen::VectorXd estimation = estim.estim;

        // auto end = std::chrono::high_resolution_clock::now();
        // if (Options::LOG_TIMER) {
        //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //     std::cout << "Estimated Points execution time: " << duration.count() << " milliseconds" << std::endl;
        // }
        // std::cout<<estimation<<std::endl;

        if(estimation.size()==0|| estimation(0)==std::numeric_limits<double>::quiet_NaN()){
            std::cout<<"Estimation point is empty"<<std::endl;
            // throw std::exception();
        }
        else{
            assert(estimation.size()==4);
            camera_variables.points.col(idx) = estimation;
        }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    if(duration.count()>1){
        // std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    }
    }

        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "Refinement point reestimation time: " << duration.count() << " Milliseconds" << std::endl;
}


void Refinement::reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras, std::vector<std::vector<int>> fixed){


        auto start= std::chrono::high_resolution_clock::now();
    for (int i =0; i<idx_cameras.size(); i++){
        int idx = (-pathway(idx_cameras(i)))-1;
        std::vector<int> fixed_points;
        fixed_points= fixed[idx_cameras[i]];
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

        // auto start = std::chrono::high_resolution_clock::now();
        EstimatedViews estim = EstimatedViews(data, camera_variables.points, fixed_points_eigen, idx , fixed_points.size());
        // auto end = std::chrono::high_resolution_clock::now();
        Eigen::VectorXd estimation= estim.estim;

        if(estimation.size()==0||estimation(0)==std::numeric_limits<double>::quiet_NaN()){
            std::cout<<"Estimation view is empty"<<std::endl;
        }
        else{
            assert(estimation.size()==12);
            Eigen::MatrixXd mat = estimation;
            mat.resize(4,3);
            camera_variables.cameras.block(idx*3,0,3,4) = mat.transpose();
        }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    if(duration.count()>1){
        // std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    }
    }
    // }

        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "Refinement view reestimation time: " << duration.count() << " Milliseconds" << std::endl;
}

