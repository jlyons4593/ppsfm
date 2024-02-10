#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma once
namespace DataStructures {

	struct InputMatrices {
	    Eigen::SparseMatrix<double> measurements;
	    Eigen::MatrixXd image_size;
	    Eigen::MatrixXd centers;
	};

	struct SfMData {
	    Eigen::MatrixXi visible; // Assuming FxN size, binary visibility matrix
	    Eigen::MatrixXd data; // Data matrix for computations
	    Eigen::MatrixXd norm_meas; // Normalized homogeneous image projection coordinates, size 3FxN
	    Eigen::MatrixXd normalisations; // Normalisation transformations, size 3Fx3
	    Eigen::MatrixXd img_meas; // Unnormalized homogeneous measurements, size 3FxN
	    Eigen::MatrixXd pinv_meas; // Pseudo-inverse or cross-product matrix, size Fx3N
	    Eigen::RowVectorXd ignored_pts; // Binary mask for ignored points, size 1xN

	    // Constructor to initialize the matrices with sizes if known at creation time
	    SfMData(int F, int N) : visible(Eigen::MatrixXi::Zero(F, N)), 
				    data(Eigen::MatrixXd::Zero(F, N)), // Modify accordingly if different dimensions are required
				    norm_meas(Eigen::MatrixXd::Zero(3*F, N)), 
				    normalisations(Eigen::MatrixXd::Zero(3*F, 3)), 
				    img_meas(Eigen::MatrixXd::Zero(3*F, N)), 
				    pinv_meas(Eigen::MatrixXd::Zero(F, 3*N)), 
				    ignored_pts(Eigen::RowVectorXd::Zero(N)) {}
	};

	struct Model {
	    Eigen::MatrixXi inliers; // FxN inliers matrix
	    Eigen::MatrixXd cameras; // 3Fx4 projective camera estimations
	    std::vector<double> timings; // Timings for reconstruction steps
	    Eigen::MatrixXd points; // 4xN projective point estimations
	    std::vector<int> pathway; // Order of views and points addition
	    std::vector<std::vector<int>> fixed; // Constraints used for views and points

	    Model(int F, int N, int k) : inliers(Eigen::MatrixXi::Zero(F, N)), 
					 cameras(Eigen::MatrixXd::Zero(3*F, 4)), 
					 points(Eigen::MatrixXd::Zero(4, N)) {
		// Initialize vectors with default sizes/values as needed
		timings.resize(4, 0.0); // Assuming there are 4 timing steps
		pathway.resize(k, 0); // Adjust 'k' based on your process
		fixed.resize(F + N); // Adjust based on your F and N
	    }
  };

	struct SfMModelSeries{
	    std::vector<Model> models;

	    // Add other fields from `data` structure as needed
	    // DataStructure data; // Assuming you have a DataStructure defined elsewhere

	    SfMModelSeries(int numModels, int F, int N, int k) {
		models.reserve(numModels);
		for (int i = 0; i < numModels; ++i) {
		    models.emplace_back(F, N, k);
		}
	    }
	};

}

