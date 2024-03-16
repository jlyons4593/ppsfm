#include "Refinement.h"
#include "dataStructures.hpp"

Refinement::Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables , int init_refine, int last_path, bool start_cameras, int type): data(data), camera_variables(camera_variables){

        
    std::cout<<camera_variables.pathway<<std::endl;
    // Eigen::VectorXi pathway_segment = camera_variables.pathway.segment(init_refine, last_path - init_refine);
}
