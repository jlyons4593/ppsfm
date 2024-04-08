#include "dataStructures.hpp"
#include "factorCompletion.h"
#include "logger.h"
#include "options.h"
#include "pairAffinityCalculator.h"
#include "ppsfm_initialisation.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

#pragma once
class PipelineManager {
private:
  /*
   * Class Variables
   *
   */
  // Variable for storing alot of the output data like visibility matrix, data
  // matrix, normalised measurements, normalisation transformation for each
  // camera stacked vertically, un normalised measurements, etc

  Eigen::SparseMatrix<double> measurements;

  // EigenMatrixx for storing the image sizes
  Eigen::MatrixXd image_size;

  DataStructures::SfMData data;
  DataStructures::ComputedCameraPoints camera_variables;
  DataStructures::ViewpairAffinity pair_affinity;

  // This function sets the many of the data variables that come from the
  // measurements passed in
  void cleanData();

  // This function sets the view_pairs and affinity measurements by making a
  // call to our measurements and affinity Generator
  void pairsAffinity();

public:
  Eigen::MatrixXd centers;
  // Constructor for PipelineManager
  PipelineManager(DataStructures::InputMatrices input);

  void runPipeline();

  // Destructor
  ~PipelineManager();
  // Function to trigger the pipeline
};
