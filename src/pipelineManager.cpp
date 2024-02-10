#include <iostream>
#include "pipelineManager.h"

PipelineManager::PipelineManager(DataStructures::InputMatrices input) 
{
    std::cout<<"Initialising Pipeline Variables"<<std::endl;
    this->measurements = input.measurements;
    this->image_size = input.image_size;
    this->centers = input.centers;

}
PipelineManager::~PipelineManager(){}
//
void PipelineManager::runPipeline(){
    std::cout<<"running pipeline"<<std::endl;
    // Prep Data
    //
    //
    // Calc pair affinity
}

void PipelineManager::cleanData(){

}
