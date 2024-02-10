#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "DataReader.cpp"
#include "dataStructures.hpp"
#include "pipelineManager.h"
#include "pipelineManager.cpp"

int main(){
    DataStructures::InputMatrices m = matread("../../../Data/sfm_dataset/bluebear.mat");
    PipelineManager manager(m);
    manager.runPipeline();
    return 0;
}
