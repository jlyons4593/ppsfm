#include "ppsfm_initialisation.h"

void Initialisation::process() {

        Logger::logSection("Starting PPSFM Initialisation");

        Logger::logSubsection("Creating Matrices for Computation");

        Eigen::MatrixXi estimated_pairs = Eigen::MatrixXi::Zero(view_pairs.rows(),view_pairs.cols());
        
        Eigen::MatrixXi not_estimated_pairs = (estimated_pairs.array() == 0).cast<int>();
        Eigen::VectorXi both_unestimated(not_estimated_pairs.rows());

        Eigen::MatrixXi visibleLogical = visible.cast<int>();

        for(int i=0; i<not_estimated_pairs.rows(); i++){
            if(not_estimated_pairs(i,0) && not_estimated_pairs(i,1)){
                both_unestimated(i) = 1;
            }
        }

        Eigen::VectorXi one_estimated(estimated_pairs.rows());

        for(int i = 0; i < estimated_pairs.rows(); i++){
            bool any_estimated = false; // Tracks if any element in the row is non-zero (estimated)
            bool all_estimated = true;  // Tracks if all elements in the row are non-zero (estimated)

            for(int j = 0; j < estimated_pairs.cols(); j++){
                if(estimated_pairs(i, j) == 0){
                    all_estimated = false; // If any element is zero, not all are estimated
                } else {
                    any_estimated = true; // If any element is non-zero, at least one is estimated
                }
            }

            // A row meets the condition if at least one element is estimated and not all are estimated
            one_estimated(i) = any_estimated && !all_estimated;
        }

        std::vector<int> sorted_idx;

        // Find indices of non-zero elements in both_unestimated
        for(int i = 0; i < both_unestimated.size(); i++) {
            if(both_unestimated(i) != 0) {
                sorted_idx.push_back(i);
            }
        }

        // Find indices of non-zero elements in one_estimated and append
        for(int i = 0; i < one_estimated.size(); i++) {
            if(one_estimated(i) != 0) {
                sorted_idx.push_back(i);
            }
        }

        // MAIN WORK LOOP

        Logger::logSubsection("Main PPSFM Initialisation Loop");

        for(int i=0; i<sorted_idx.size(); i++){
        // for(int i=0; i<1; i++){
            int first_view = view_pairs(sorted_idx[i], 0); // Remember, C++ uses 0-based indexing
            int second_view = view_pairs(sorted_idx[i], 1);

            // logical and for creating visible points rowVector
            auto logicalAnd = [](int a, int b) -> int {
                return a && b;
            };
            Eigen::VectorXi visible_points = visible.row(first_view).binaryExpr(visible.row(second_view), logicalAnd);
            // std::cout<<visible_points.size()<<std::endl;
            std::vector<int> first_view_idx;
            std::vector<int> second_view_idx;
            for(int i = 3*first_view; i < 3*first_view+2; i++) {
                first_view_idx.push_back(i);
            }
            for(int i = 3*second_view; i < 3*second_view+2; i++) {
                second_view_idx.push_back(i);
            }
            int number_of_visible = (visible_points.array() != 0).count();

            Eigen::MatrixXd block1(second_view_idx.size(),number_of_visible);
            Eigen::MatrixXd block2(first_view_idx.size(),number_of_visible);
            
            //Fill measurement block 1
            for (int i = 0; i < second_view_idx.size(); ++i) {
                for (int j = 0; j < visible_points.size(); ++j) {
                    if(visible_points(j)>0){
                        block1(i, j) = measurement(second_view_idx[i], j);
                    }
                }
            }
            // Fill measuremnent block 2
            for (int i = 0; i < first_view_idx.size(); ++i) {
                for (int j = 0; j < visible_points.size(); ++j) {
                    if(visible_points(j)>0){
                        block2(i, j) = measurement(first_view_idx[i], j);
                    }
                }
            }

            // Call to estimated fundamental matrix
            EstimatedFundamentalMatrix estimated_fundamental_matrix(block1,block2, 99.99, 1000, 1e-3);

            Eigen::Array<bool,Eigen::Dynamic,1> inliers = estimated_fundamental_matrix.getInliers();
            Eigen::MatrixXd fundamental_matrix= estimated_fundamental_matrix.getFundamentalMatrix();

            int status=0;

            if (status==0){
                std::vector<int> initial_views = {first_view, second_view};

                std::vector<int> visible_idx;
                for(int i = 0; i < visible_points.size(); ++i) {
                    if(visible_points(i)) {
                        visible_idx.push_back(i); // Store zero-based indices of visible points
                    }
                }

                std::vector<int> initial_points;
                for(int i = 0; i < inliers.size(); ++i) {
                    if(inliers(i)) { // If the inlier_idx at position i is true, select the corresponding index
                        initial_points.push_back(visible_idx[i]);
                    }
                }
                camera_variables = compute_cams(initial_views, initial_points, fundamental_matrix);
                break;
            }
        }
        Logger::logSection("End of Camera Var Computation");

    }
