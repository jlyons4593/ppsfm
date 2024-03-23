#include "Refinement.h"
#include "dataStructures.hpp"
#include "options.h"
#include "pyramidalVisibilityScore.h"
#include <variant>

Refinement::Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables , int init_refine, int last_path, bool start_cameras, int type): data(data), camera_variables(camera_variables){

        
    Eigen::VectorXi pathway_segment = camera_variables.pathway.segment(init_refine-1, last_path+1 );

    Eigen::Array<std::variant<int, std::vector<int>, Eigen::VectorXi>, Eigen::Dynamic, 1> fixed = camera_variables.fixed.segment(init_refine-1, last_path+1);

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
    Eigen::VectorXi idx_cameras(camera_variables.negative_pathway.size()); 

    int camera_count =0;
    int point_count =0;
    
    Eigen::VectorXi known_points = pathway_segment(idx_points);

    for(int i =0; i<pathway_segment.size();i++){
        if(pathway_segment(i)<0){
            idx_cameras(camera_count) = i;
            camera_count++;
        }else{
            idx_points(point_count) = i;
            point_count++;
        }

    }
    // Initialize idx_cameras with appropriate values


    Eigen::VectorXi selected_pathway = pathway_segment(idx_cameras);
    Eigen::VectorXi negated_scaled_pathway = selected_pathway * -3;

    Eigen::MatrixXi replicated_pathway = negated_scaled_pathway.replicate(1, 3);
    replicated_pathway.resizeLike(Eigen::MatrixXi::Ones(idx_cameras.size(), 3));

    Eigen::RowVectorXi sequence(3);
    sequence << 2, 1, 0;

    Eigen::MatrixXi replicated_sequence = sequence.replicate(idx_cameras.size(), 1);
    replicated_pathway = replicated_pathway.array()-1;


    Eigen::VectorXi known_cams3 = (replicated_pathway - replicated_sequence).reshaped();


    std::vector<int> rowSize;
    for (int i =0; i<data.visible.cols(); i++){
        rowSize.push_back(i);
    }
    Eigen::VectorXi positive_path = camera_variables.positive_pathway;

    std::sort(rowSize.begin(), rowSize.end());
    std::sort(positive_path.begin(), positive_path.end());

    // // Compute set difference
    std::vector<int> diff;
    std::set_difference(rowSize.begin(), rowSize.end(), positive_path.begin(), positive_path.end(), std::back_inserter(diff));
        

    for(int i=0; i<diff.size(); i++){
        data.visible.col(diff[i]).array() = 0;
    }

    std::vector<int> colSize;
    for (int i =0; i<data.visible.rows(); i++){
        colSize.push_back(i); 
    }
    Eigen::VectorXi negative_path = camera_variables.negative_pathway;
    std::vector<int> diff2;
    std::set_difference(colSize.begin(), colSize.end(), negative_path.begin(), negative_path.end(), std::back_inserter(diff2));

    for(int i=0; i<diff2.size(); i++){
        data.visible.row(diff2[i]).array() = 0;
    }

    int iteration_number = 0 ;
    double change_camera = std::numeric_limits<double>::infinity();

    double change_points = std::numeric_limits<double>::infinity();


    // while(iteration_number<max_iter){
    while(iteration_number<1){
        Eigen::MatrixXd old_cameras(known_cams3.size(), camera_variables.cameras.cols());
        for (int i = 0; i < known_cams3.size(); ++i) {
            old_cameras.row(i) = camera_variables.cameras.row(known_cams3(i));

        }

        Eigen::MatrixXd old_points = camera_variables.points(Eigen::all, known_points);

        // if(start_cameras){
            // camera_variables.cameras = reestimate_all_views(pathway_segment, idx_cameras);
            // camera_variables.points = reestimate_all_views(pathway_segment, idx_points);

        // }
        // else{
            // camera_variables.points = reestimate_all_points(pathway_segment, idx_points);
        reestimate_all_views(pathway_segment, idx_cameras);
        // }
        iteration_number++;
    }


}

    void Refinement::reestimate_all_views(Eigen::VectorXi pathway, Eigen::VectorXi idx_cameras){

        // for(int j =0; j<idx_cameras.size(); j++){
        for(int j =0; j<1; j++){
            int idx = (-pathway(idx_cameras(j)))-1;
            std::cout<<idx<<std::endl;
            std::variant<int, std::vector<int>, Eigen::VectorXi> fixed_points = fixed[idx_cameras(j)];
            Eigen::VectorXd visible_points= data.visible.row(idx);
            

            // Access the value based on its type
            // if (std::holds_alternative<int>(fixed_points)) {
                // int fixed_views_value = std::get<int>(fixed_points);
                // Use fixed_views_value as needed
            // } 
            // else if (std::holds_alternative<std::vector<int>>(fixed_points)) {
            //     std::vector<int> fixed_views_vec = std::get<std::vector<int>>(fixed_points);
            //     // Use fixed_views_vec as needed
            // } else if (std::holds_alternative<Eigen::VectorXi>(fixed_points)) {
            //     Eigen::VectorXi fixed_views_eigen = std::get<Eigen::VectorXi>(fixed_points);
            //     // Use fixed_views_eigen as needed
            // }

            // Please find a cleaner way to deal with this variant

            // if(fixed_std_vec->size()!=0){
            //     std::cout<<fixed_std_vec<<std::endl;
            //     for(int i=0; i<fixed_std_vec->size(); i++ ){
            //         fixed_std_vec->at(i) = 0;
            //     }
            // }
            // else if(fixed_eig_vec->size()!=0){
            //     std::cout<<fixed_eig_vec<<std::endl;
            //     for(int i=0; i<fixed_eig_vec->size(); i++ ){
            //         (*fixed_eig_vec)(i) = 0;
            //     }
            // }
            // else{
            //     std::cout<<*fixed_int<<std::endl;
            //     (*fixed_int) = 0;
            // }
            //

        }


    }
    Eigen::MatrixXd Refinement::reestimate_all_points(Eigen::VectorXi pathway, Eigen::VectorXi idx_points){ 

    }
