#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include "DataReader.cpp"
#include "dataStructures.hpp"
#include "pipelineManager.h"
#include "pipelineManager.cpp"
    
int main(){
    DataStructures::InputMatrices matrices = matread("../../../Data/sfm_dataset/bluebear.mat");
    PipelineManager manager(matrices);
    manager.runPipeline();
    return 0;
}
