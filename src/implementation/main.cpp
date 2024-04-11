#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "DataReader.hpp"
#include "dataStructures.hpp"
#include "pipelineManager.h"
#include <omp.h>
 
int main(){
    DataStructures::InputMatrices matrices = DataReader::matread("../../../Data/sfm_dataset/dino4983.mat");
    // DataStructures::InputMatrices matrices = DataReader::matread("../../../Data/sfm_dataset/cherub2Colmap.mat");
    PipelineManager manager(matrices);
    manager.runPipeline();
    return 0;
}
