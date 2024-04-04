#include "DataWriter.hpp"
#include <iostream>
#include "pipelineManager.h"
#include "Refinement.h"
#include "dataCleaning.h"
#include "dataStructures.hpp"

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

    Logger::logSection("Finding Initial Sub-Problem");

    std::unique_ptr<Initialisation> initialiser = std::make_unique<Initialisation>(data.normalised_measurements,data.visible, pair_affinity.view_pairs, pair_affinity.Affinity);
    initialiser->process();
    camera_variables = initialiser->getCameraPoints();
    Eigen::MatrixXd initial_cameras= camera_variables.cameras;
    Eigen::MatrixXd initial_points= camera_variables.points;
    Eigen::VectorXi initial_pathway=camera_variables.pathway;
    std::vector<std::vector<int>> initial_fixed = camera_variables.fixed;

  for(std::vector<int>& i: initial_fixed){
      for(int& j: i){
          j += 1;
      }

  }
    std::unique_ptr<FactorCompletion> completion = std::make_unique<FactorCompletion>(data, camera_variables,pair_affinity, image_size, centers);
    completion->process();
    camera_variables = completion->getCameraVariables();
    data = completion->getData();
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> inliers = completion->getInliers();

    if(Options::FINAL_REFINEMENT){
        Refinement refinement = Refinement(data, camera_variables, inliers, camera_variables.pathway, camera_variables.fixed, false, 2);
        camera_variables.cameras = refinement.getCameras();
        camera_variables.points= refinement.getPoints();
    }
     
    std::vector<std::vector<int>> fixed = camera_variables.fixed;
  for(std::vector<int>& i: fixed){
      for(int& j: i){
          j += 1;
      }

  }
    std::vector<int> positiveValues;
    std::copy_if(
            camera_variables.pathway.data(),
            camera_variables.pathway.data() + camera_variables.pathway.size(),
            std::back_inserter(positiveValues), [](int value) { return value > 0; });
    // Now create a new Eigen vector from the std::vector of positive values
    Eigen::Map<Eigen::VectorXi> valid_points(positiveValues.data(),
            positiveValues.size());
    std::vector<int> negativeValues;
    std::copy_if(
            camera_variables.pathway.data(),
            camera_variables.pathway.data() + camera_variables.pathway.size(),
            std::back_inserter(negativeValues), [](int value) { return value < 0; });
    // Now create a new Eigen vector from the std::vector of positive values
    Eigen::Map<Eigen::VectorXi> valid_cams(negativeValues.data(),
            negativeValues.size());

    Eigen::VectorXi estimated_points = valid_points.array() - 1;
    Eigen::VectorXi estimated_views= -(valid_cams.array()) - 1;

    std::vector<double> timings;

    DataStructures::ModelData model;


    model.initial_cameras=initial_cameras;
    model.initial_points=initial_points;
    model.initial_fixed=initial_fixed;
    model.initial_pathway=initial_pathway;
    model.cameras=camera_variables.cameras;
    model.points=camera_variables.points;
    model.pathway=camera_variables.pathway;
    model.fixed=fixed;
    model.inliers= inliers;
    
    DataStructures::SfMData finalData;
    finalData.cost_function_data = data.cost_function_data;
    finalData.normalisations= data.normalisations;
    finalData.normalised_measurements = data.normalised_measurements;
    finalData.visible = data.visible;
    finalData.image_measurements= data.image_measurements;
    finalData.pseudo_inverse_measurements= data.pseudo_inverse_measurements;
    finalData.removed_points= data.removed_points;
    std::cout<<"Data has been saved to its structures"<<std::endl;

    bool save =true;
    if(save){
    WriteStructsToMatFile(model, data, "bluebear.mat");
    std::cout<<"Data has been saved to its matfile"<<std::endl;
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
