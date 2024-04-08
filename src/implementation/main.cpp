#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "DataReader.hpp"
#include "dataStructures.hpp"
#include "pipelineManager.h"
    
int main(){
    DataStructures::InputMatrices matrices = DataReader::matread("../../../Data/sfm_dataset/dino4983.mat");

    PipelineManager manager(matrices);
    manager.runPipeline();
    return 0;
}
