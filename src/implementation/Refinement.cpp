#include "Refinement.h"
#include "dataStructures.hpp"
#include "options.h"

Refinement::Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables , int init_refine, int last_path, bool start_cameras, int type): data(data), camera_variables(camera_variables){

        
    Eigen::VectorXi pathway_segment = camera_variables.pathway.segment(init_refine-1, last_path+1 );

    Eigen::Array<std::variant<int, std::vector<int>, Eigen::VectorXi>, Eigen::Dynamic, 1> fixed_segment = camera_variables.fixed.segment(init_refine-1, last_path+1);

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
        // std::cout<<point_count<<std::endl;
        // std::cout<<camera_count<<std::endl;
        // std::cout<<i<<std::endl;

    }
    std::cout<<idx_cameras<<std::endl;
// Initialize idx_cameras with appropriate values

    Eigen::VectorXi selected_pathway = pathway_segment(idx_cameras);
    Eigen::VectorXi negated_scaled_pathway = selected_pathway * 3;

    Eigen::MatrixXi replicated_pathway = negated_scaled_pathway.replicate(1, 3);
    replicated_pathway.resizeLike(Eigen::MatrixXi::Ones(idx_cameras.size(), 3));

    Eigen::RowVectorXi sequence(3);
    sequence << 2, 1, 0;

    Eigen::MatrixXi replicated_sequence = sequence.replicate(idx_cameras.size(), 1);

    Eigen::VectorXi known_cams3 = (replicated_pathway - replicated_sequence).reshaped();
    

    // Eigen::VectorXi known_points = pathway_segment(idx_points);
    // std::cout<<idx_cameras<<std::endl;
}
