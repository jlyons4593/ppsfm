
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

void PipelineManager::runPipeline(){
    cleanData();    
    pairsAffinity();
    for(int i =0; i<Options::MAX_MODELS; i++){

        Logger::logSection("Finding Initial Sub-Problem");

        std::unique_ptr<Initialisation> initialiser = std::make_unique<Initialisation>(data.normalised_measurements,data.visible, pair_affinity.view_pairs, pair_affinity.Affinity);
        initialiser->process();
        camera_variables = initialiser->getCameraPoints();

        std::unique_ptr<FactorCompletion> completion = std::make_unique<FactorCompletion>(data, camera_variables,pair_affinity, image_size, centers);
        completion->process();



    }
}



void PipelineManager::cleanData(){
    std::unique_ptr<DataCleaningStage> dataCleaner(new DataCleaningStage); 
    dataCleaner->process(measurements);
    this->data = dataCleaner->getData();
}
void PipelineManager::pairsAffinity(){
    std::unique_ptr<PairAffinityCalculator> pairHandler(new PairAffinityCalculator); 
    pairHandler->process(data.image_measurements, data.visible, image_size);
    this->pair_affinity = pairHandler->getPairAffinity();

}
