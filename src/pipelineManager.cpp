#include <iostream>
#include "pipelineManager.h"
#include "dataCleaning.h"

PipelineManager::PipelineManager(DataStructures::InputMatrices input) 
{
    std::cout<<"Initialising Pipeline Variables"<<std::endl;
    this->measurements = input.measurements;
    this->image_size = input.image_size;
    this->centers = input.centers;

}
PipelineManager::~PipelineManager(){}
//
// void PipelineManager::runPipeline(){
    // std::cout<<"running pipeline"<<std::endl;
    // Prep Data
    //
    //
    // Calc pair affinity
// }

void PipelineManager::cleanData(){
    std::unique_ptr<DataCleaningStage> dataCleaner(new DataCleaningStage); 
    dataCleaner->process(measurements);
    this->data = dataCleaner->getData();
}
void PipelineManager::pairsAffinity(){
    std::unique_ptr<PairAffinityStage> pairHandler(new PairAffinityStage); 
    pairHandler->process(data.image_measurements, data.visible, image_size);
    this->pair_affinity = pairHandler->getPairAffinity();

}

