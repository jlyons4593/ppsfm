#include <Eigen/Dense>
#include "ppsfm_initialisation.h"
#include "logger.h"
#include "options.h"
#include <functional>
#include <iostream>
#include <Eigen/Sparse>
#include "dataStructures.hpp"
#include "pairAffinityCalculator.h"
    
// This class takes in measurements, img_sizes, centers, options
//     * measurements: Original image coordinates (2FxN sparse matrix where missing data are [0;0])
//     * image_size: Size of the images, used to compute the pyramidal visibility scores (2xF)
//     * centers: Principal points coordinates for each camera (2xF)
//     * options:  Structure containing options (must be initialized by ppsfm_options to contains all necessary fields)
// It will then manage the models and data and be able to return these variables once properly instantiated.
//   models: A cell of structures (one for each model) with the following fields:
//         - inliers: Inliers matrix binary mask (FxN)
//         - cameras: Projective estimation of the cameras that has been completed (3Fx4)
//         - timings: Timings of the four steps of the reconstruction
//         - points: Projective estimation of the points that has been completed (4xN)
//         - pathway: Array containing the order in which views (negative) and points (positive) has been added (1 x k, k <= F+N)
//         - fixed: Cell containing arrays of the points or views used in the constraints when adding views and points (1 x F+N) 
// * data: Structure containing all data in the following fields:
//         - visible: Visibility matrix binary mask (FxN)
//         - data: Data matrix used for computations
//         - norm_meas: Normalized homogeneous image projection coordinates (3FxN)
//         - normalisations: Normalisation transformation for each camera stacked vertically (3Fx3)
//         - img_meas: Unnormaliazed homogeneous measurements of the projections, used for computing errors and scores (3FxN)
//         - pinv_meas: Matrix containing the data for elimination of projective depths, cross-product matrix or pseudo-inverse of the of the normalized homogeneous coordinates (Fx3N)
//         - ignored_pts: Binary mask indicated points removed from the reconstruction due to not having enought projections (1xN)

#pragma once
class PipelineManager{
private:
    /*
     * Class Variables
     * 
     */
    // Variable for storing alot of the output data like visibility matrix, data matrix, normalised measurements, normalisation transformation for each camera stacked vertically, un normalised measurements, etc



    Eigen::SparseMatrix<double> measurements;
        
    // EigenMatrixx for storing the image sizes
    Eigen::MatrixXd image_size;

    Eigen::MatrixXd centers;
    DataStructures::SfMData data;
    DataStructures::ComputedCameraPoints camera_variables;
    DataStructures::ViewpairAffinity pair_affinity;

    // This function sets the many of the data variables that come from the measurements passed in
    void cleanData();
     
    // This function sets the view_pairs and affinity measurements by making a call to our measurements and affinity Generator
    void pairsAffinity();

    
public:

    // Constructor for PipelineManager
    PipelineManager(DataStructures::InputMatrices input);

    void runPipeline() {
        cleanData();    
        pairsAffinity();
        for(int i =0; i<Options::MAX_MODELS; i++){

            Logger::logSection("Finding Initial Sub-Problem");
            
            std::unique_ptr<Initialisation> initialiser = std::make_unique<Initialisation>(data.normalised_measurements,data.visible, pair_affinity.view_pairs, pair_affinity.Affinity);
            initialiser->process();
            camera_variables = initialiser->getCameraPoints();


        }

    }

    // Destructor
    ~PipelineManager();
    // Function to trigger the pipeline

};
