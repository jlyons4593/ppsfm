#include "dataStructures.hpp"
#include "mat.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

inline void WriteStructsToMatFile(const DataStructures::ModelData& modelData, const DataStructures::SfMData& sfmData, const std::string& fileName) {
    // Create a new MAT-file for writing.
    MATFile *pmat = matOpen(fileName.c_str(), "w");
    if (pmat == NULL) {
      std::cout<<"MATLAB:writeStructsToMatFile:fileCreationFailed"<< " Unable to create MAT-file."<<std::endl;
        return;
    }
    // Close the MAT-file.
    mxArray *initial_cameras = mxCreateDoubleMatrix(modelData.initial_cameras.rows(), modelData.initial_cameras.cols(), mxREAL);
    std::copy(modelData.initial_cameras.data(), modelData.initial_cameras.data() + modelData.initial_cameras.size(), mxGetPr(initial_cameras));

    // Write the mxArray to the MAT-file.
    if (matPutVariable(pmat, "initial_cameras", initial_cameras) != 0) {
      std::cout<<"MATLAB:writeStructsToMatFile:saveFailed. Unable to save 'initial_cameras' to MAT-file."<<std::endl;;
        mxDestroyArray(initial_cameras);
        matClose(pmat);
        return;
    }

    // Continuing to save other fields in ModelData
    // Convert and save 'initial_points'
    mxArray *initial_points = mxCreateDoubleMatrix(modelData.initial_points.rows(), modelData.initial_points.cols(), mxREAL);
    std::copy(modelData.initial_points.data(), modelData.initial_points.data() + modelData.initial_points.size(), mxGetPr(initial_points));
    if (matPutVariable(pmat, "initial_points", initial_points) != 0) {
      std::cout<<"Unable to save 'initial_points' to MAT-file."<<std::endl;
        mxDestroyArray(initial_points);
        matClose(pmat);
        return;
    }

    // The process should be repeated for all fields, ensuring correct data type conversion.
    // For Eigen::VectorXi, we need to convert it to double for mxArray, as mxCreateNumericMatrix does not support int types directly.
    mxArray *initial_pathway = mxCreateDoubleMatrix(modelData.initial_pathway.size(), 1, mxREAL);
    // Need to convert int vector to double
    Eigen::VectorXd initial_pathway_double = modelData.initial_pathway.cast<double>();
    std::copy(initial_pathway_double.data(), initial_pathway_double.data() + initial_pathway_double.size(), mxGetPr(initial_pathway));
    if (matPutVariable(pmat, "initial_pathway", initial_pathway) != 0) {
      std::cout<<"Unable to save 'initial_pathway' to MAT-file."<<std::endl;
        mxDestroyArray(initial_pathway);
        matClose(pmat);
        return;
    }

    // Example for saving std::vector<std::vector<int>> requires converting to a MATLAB cell array
    // For simplicity, this example assumes 'initial_fixed' contains data. Add checks in real scenarios.
    mxArray *initial_fixed = mxCreateCellMatrix(modelData.initial_fixed.size(), 1);
    for (size_t i = 0; i < modelData.initial_fixed.size(); ++i) {
        mxArray *row = mxCreateDoubleMatrix(modelData.initial_fixed[i].size(), 1, mxREAL);
        for (size_t j = 0; j < modelData.initial_fixed[i].size(); ++j) {
            // Convert int to double for mxArray
            mxGetPr(row)[j] = static_cast<double>(modelData.initial_fixed[i][j]);
        }
        mxSetCell(initial_fixed, i, row);
    }
    if (matPutVariable(pmat, "initial_fixed", initial_fixed) != 0) {
      std::cout<<"Unable to save 'initial_fixed' to MAT-file."<<std::endl;
        mxDestroyArray(initial_fixed);
        matClose(pmat);
        return;
    }



    mxArray *cameras = mxCreateDoubleMatrix(modelData.cameras.rows(), modelData.cameras.cols(), mxREAL);
    std::copy(modelData.cameras.data(), modelData.cameras.data() + modelData.cameras.size(), mxGetPr(cameras));
    if (matPutVariable(pmat, "cameras", cameras) != 0) {
      std::cout<<"Unable to save 'cameras' to MAT-file."<<std::endl;
        mxDestroyArray(cameras);
        matClose(pmat);
        return;
    }

    // Convert and save 'points'
    mxArray *points = mxCreateDoubleMatrix(modelData.points.rows(), modelData.points.cols(), mxREAL);
    std::copy(modelData.points.data(), modelData.points.data() + modelData.points.size(), mxGetPr(points));
    if (matPutVariable(pmat, "points", points) != 0) {
      std::cout<<"Unable to save 'points' to MAT-file."<<std::endl;
        mxDestroyArray(points);
        matClose(pmat);
        return;
    }

    // Convert and save 'pathway'
    mxArray *pathway = mxCreateDoubleMatrix(modelData.pathway.size(), 1, mxREAL);
    Eigen::VectorXd pathway_double = modelData.pathway.cast<double>();
    std::copy(pathway_double.data(), pathway_double.data() + pathway_double.size(), mxGetPr(pathway));
    if (matPutVariable(pmat, "pathway", pathway) != 0) {
      std::cout<<"Unable to save 'pathway' to MAT-file."<<std::endl;
        mxDestroyArray(pathway);
        matClose(pmat);
        return;
    }

    // Convert and save 'fixed' similar to 'initial_fixed'
    mxArray *fixed = mxCreateCellMatrix(modelData.fixed.size(), 1);
    for (size_t i = 0; i < modelData.fixed.size(); ++i) {
        mxArray *row = mxCreateDoubleMatrix(modelData.fixed[i].size(), 1, mxREAL);
        for (size_t j = 0; j < modelData.fixed[i].size(); ++j) {
            mxGetPr(row)[j] = static_cast<double>(modelData.fixed[i][j]);
        }
        mxSetCell(fixed, i, row);
    }
    if (matPutVariable(pmat, "fixed", fixed) != 0) {
      std::cout<<"Unable to save 'fixed' to MAT-file."<<std::endl;
        mxDestroyArray(fixed);
        matClose(pmat);
        return;
    }

    // Convert and save 'inliers'
    mxArray *inliers = mxCreateLogicalMatrix(modelData.inliers.rows(), modelData.inliers.cols());
    memcpy(mxGetLogicals(inliers), modelData.inliers.data(), modelData.inliers.size() * sizeof(mxLogical));
    if (matPutVariable(pmat, "inliers", inliers) != 0) {
      std::cout<<"Unable to save 'inliers' to MAT-file."<<std::endl;
        mxDestroyArray(inliers);
        matClose(pmat);
        return;
    }

    // Now, converting and saving SfMData
    // This follows the same pattern as above for each field
    // Below is an example for 'visible'. Repeat for all fields in SfMData.
    mxArray *visible = mxCreateDoubleMatrix(sfmData.visible.rows(), sfmData.visible.cols(), mxREAL);
    std::copy(sfmData.visible.data(), sfmData.visible.data() + sfmData.visible.size(), mxGetPr(visible));
    if (matPutVariable(pmat, "visible", visible) != 0) {
      std::cout<<"Unable to save 'visible' to MAT-file."<<std::endl;
        mxDestroyArray(visible);
        matClose(pmat);
        return;
    }
   // Convert and save 'cost_function_data'
    mxArray *cost_function_data = mxCreateDoubleMatrix(sfmData.cost_function_data.rows(), sfmData.cost_function_data.cols(), mxREAL);
    std::copy(sfmData.cost_function_data.data(), sfmData.cost_function_data.data() + sfmData.cost_function_data.size(), mxGetPr(cost_function_data));
    if (matPutVariable(pmat, "cost_function_data", cost_function_data) != 0) {
      std::cout<<"Unable to save 'cost_function_data' to MAT-file."<<std::endl;
        mxDestroyArray(cost_function_data);
        matClose(pmat);
        return;
    }

    // Convert and save 'normalised_measurements'
    mxArray *normalised_measurements = mxCreateDoubleMatrix(sfmData.normalised_measurements.rows(), sfmData.normalised_measurements.cols(), mxREAL);
    std::copy(sfmData.normalised_measurements.data(), sfmData.normalised_measurements.data() + sfmData.normalised_measurements.size(), mxGetPr(normalised_measurements));
    if (matPutVariable(pmat, "normalised_measurements", normalised_measurements) != 0) {
      std::cout<<"Unable to save 'normalised_measurements' to MAT-file."<<std::endl;
        mxDestroyArray(normalised_measurements);
        matClose(pmat);
        return;
    }

    // Convert and save 'normalisations'
    mxArray *normalisations = mxCreateDoubleMatrix(sfmData.normalisations.rows(), sfmData.normalisations.cols(), mxREAL);
    std::copy(sfmData.normalisations.data(), sfmData.normalisations.data() + sfmData.normalisations.size(), mxGetPr(normalisations));
    if (matPutVariable(pmat, "normalisations", normalisations) != 0) {
      std::cout<<"Unable to save 'normalisations' to MAT-file."<<std::endl;
        mxDestroyArray(normalisations);
        matClose(pmat);
        return;
    }

    // Convert and save 'image_measurements'
    mxArray *image_measurements = mxCreateDoubleMatrix(sfmData.image_measurements.rows(), sfmData.image_measurements.cols(), mxREAL);
    std::copy(sfmData.image_measurements.data(), sfmData.image_measurements.data() + sfmData.image_measurements.size(), mxGetPr(image_measurements));
    if (matPutVariable(pmat, "image_measurements", image_measurements) != 0) {
      std::cout<<"Unable to save 'image_measurements' to MAT-file."<<std::endl;
        mxDestroyArray(image_measurements);
        matClose(pmat);
        return;
    }

    // Convert and save 'pseudo_inverse_measurements'
    mxArray *pseudo_inverse_measurements = mxCreateDoubleMatrix(sfmData.pseudo_inverse_measurements.rows(), sfmData.pseudo_inverse_measurements.cols(), mxREAL);
    std::copy(sfmData.pseudo_inverse_measurements.data(), sfmData.pseudo_inverse_measurements.data() + sfmData.pseudo_inverse_measurements.size(), mxGetPr(pseudo_inverse_measurements));
    if (matPutVariable(pmat, "pseudo_inverse_measurements", pseudo_inverse_measurements) != 0) {
      std::cout<<"Unable to save 'pseudo_inverse_measurements' to MAT-file."<<std::endl;
        mxDestroyArray(pseudo_inverse_measurements);
        matClose(pmat);
        return;
    }

    // Convert and save 'removed_points' - note this is a logical array
    mxArray *removed_points = mxCreateLogicalMatrix(1, sfmData.removed_points.size());
    Eigen::Array<bool, 1, Eigen::Dynamic>::Map(mxGetLogicals(removed_points), sfmData.removed_points.size()) = sfmData.removed_points;
    if (matPutVariable(pmat, "removed_points", removed_points) != 0) {
      std::cout<<"Unable to save 'removed_points' to MAT-file."<<std::endl;
        mxDestroyArray(removed_points);
        matClose(pmat);
        return;
    }

    // After saving all fields, remember to clean up by destroying the created mxArrays
    mxDestroyArray(initial_points);
    mxDestroyArray(initial_cameras);
    mxDestroyArray(initial_pathway);
    mxDestroyArray(initial_fixed);
    mxDestroyArray(cameras);
    mxDestroyArray(points);
    mxDestroyArray(pathway);
    mxDestroyArray(fixed);
    mxDestroyArray(inliers);

    mxDestroyArray(visible);
    mxDestroyArray(cost_function_data);
    mxDestroyArray(normalised_measurements);
    mxDestroyArray(normalisations);
    mxDestroyArray(image_measurements);
    mxDestroyArray(pseudo_inverse_measurements);
    mxDestroyArray(removed_points);

    if (matClose(pmat) != 0) {
      std::cout<<"MATLAB:writeStructsToMatFile:fileCloseFailed. Unable to close MAT-file."<<std::endl;
        return;
    }
}
