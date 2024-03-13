#include "estimatedFundamentalMatrix.h"

Eigen::VectorXd EstimatedFundamentalMatrix::linvec(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    Eigen::VectorXd v(9); // Initialize a vector of size 9

    v(0) = b(0) * a(0);
    v(1) = b(0) * a(1);
    v(2) = b(0);
    v(3) = b(1) * a(0);
    v(4) = b(1) * a(1);
    v(5) = b(1);
    v(6) = a(0);
    v(7) = a(1);
    v(8) = 1;

    return v;
}

std::vector<int> EstimatedFundamentalMatrix::get_random_sequence(int number_of_projections){

    std::vector<int> test_set(number_of_projections);
    std::iota(test_set.begin(), test_set.end(), 1); // Fill with 1, 2, ..., num_projs

    // Create a random number generator

    // Shuffle the vector to get a random permutation
    std::shuffle(test_set.begin(), test_set.end(), rng);
    if (test_set.size() > 8) {
        test_set.resize(8);
    }
    return test_set;
}

    Eigen::MatrixXd EstimatedFundamentalMatrix::estimate(Eigen::MatrixXd coeffs, bool enforce_rank){
        // Using jacobi SVD as faster for smaller problems
        // Possible performance improvement if I can utilise this as a thin V instead of a full V as will use less compute but will become inconsistent with matlab
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(coeffs, Eigen::ComputeFullV);

        // Getting the V value
        Eigen::MatrixXd fun_vec = svd.matrixV();

        // Assuming the last column of fun_vec is what we want to reshape into fun_mat
        Eigen::VectorXd last_col = fun_vec.col(fun_vec.cols()-1);

        Eigen::MatrixXd fun_mat = last_col.reshaped(3,3);

        if(enforce_rank){
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(fun_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::MatrixXd U = svd.matrixU();
            Eigen::MatrixXd V = svd.matrixV();
            Eigen::VectorXd S = svd.singularValues();

            // Set the smallest singular value to 0 to enforce a rank of 2
            S(S.size() - 1) = 0;

            // Recompose the matrix with the modified singular values
            Eigen::MatrixXd S_modified = Eigen::MatrixXd::Zero(fun_mat.rows(), fun_mat.cols());
            for(int i = 0; i < S.size(); ++i) {
                S_modified(i, i) = S(i);
            }

            fun_mat = U * S_modified * V.transpose();
        }
        return fun_mat;
        // DataStructures::printColsRows(fun_mat, "fundamental_matrix");
    }

    std::pair<Eigen::MatrixXd,Eigen::Array<bool, Eigen::Dynamic, 1>> EstimatedFundamentalMatrix::estimateRansac(Eigen::MatrixXd coeffs, double confidence, int max_iter, double dist_thresh){
        
        int number_of_projections = coeffs.rows();
        double best_score = std::numeric_limits<double>::infinity();
        Eigen::MatrixXd best_estim = Eigen::MatrixXd::Zero(3,3); // Assuming best_estim should be a dynamic-size vector; adjust if needed
        Eigen::Array<bool, Eigen::Dynamic, 1> best_inliers= Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(number_of_projections); // Initializes a vector of size num_projs with all elements set to false
        int number_iterator = 0;
        int stop_iterator = max_iter;
        double log_conf = std::log(1 - confidence / 100.0);

        while(number_iterator<=stop_iterator){
            std::vector<int> test_set = get_random_sequence(number_of_projections);
            Eigen::MatrixXd coeff_block(test_set.size(), coeffs.cols());
            for (size_t i = 0; i < test_set.size(); ++i) {
                coeff_block.row(i) = coeffs.row(test_set[i]);
            }
            Eigen::MatrixXd estimation = estimate(coeff_block,false);

            Eigen::VectorXd estimationVec = Eigen::Map<Eigen::VectorXd>(estimation.data(), estimation.size());
            Eigen::VectorXd product = coeff_block * estimationVec;

            if ((product.norm() / std::sqrt(8)) < dist_thresh) {
                Eigen::VectorXd errors = (coeffs * estimationVec).array().square();
                Eigen::Array<bool, Eigen::Dynamic, 1> inliers = (errors.array() < dist_thresh);
                double score = 0.0;
                int inlier_count = 0;
                for (int i = 0; i < errors.size(); ++i) {
                    if (inliers(i)) {
                        score += errors(i); // Sum errors for inliers
                        inlier_count++; // Count inliers
                    }
                }
                score += (inliers.size() - inlier_count) * dist_thresh; 
                if(inliers.count()>8 && score<7*best_score){

                    std::vector<int> trueIndices;
                    for (int i = 0; i < inliers.size(); ++i) {
                        if (inliers(i)) {
                            trueIndices.push_back(i);
                        }
                    }

                    Eigen::MatrixXd block(trueIndices.size(), coeffs.cols());
                    for (size_t i = 0; i < trueIndices.size(); ++i) {
                        block.row(i) = coeffs.row(trueIndices[i]);
                    }
                    estimation = estimate(block, false);
                    errors = (block * estimationVec).array().square();
                    score = 0.0;
                    inlier_count=0;
                    for (int i = 0; i < errors.size(); ++i) {
                        if (inliers(i)) {
                            score += errors(i); // Sum errors for inliers
                            inlier_count++; // Count inliers
                        }
                    }
                    score += (inliers.size() - inlier_count) * dist_thresh; 
                    if(score<best_score){
                        best_score =score;
                        best_estim = estimation;
                        best_inliers = inliers;

                        double ratio = double(inlier_count)/number_of_projections;
                        double eps = std::numeric_limits<double>::epsilon();
                        double prob = std::max(eps, std::min(1.0 - eps, 1.0 - std::pow(ratio, 8)));
                        stop_iterator = std::min(static_cast<int>(std::ceil(log_conf / std::log(prob))), max_iter);
                    }
                }
            }
            number_iterator++;
        }
        return {best_estim, best_inliers};
    }

Eigen::MatrixXd EstimatedFundamentalMatrix::applyLinvecToProjections(const Eigen::MatrixXd& projs1, const Eigen::MatrixXd& projs2) {
    if (projs1.cols() != projs2.cols()) {
        throw std::invalid_argument("projs1 and projs2 must have the same number of columns.");
    }

    Eigen::MatrixXd coeffs(projs1.cols(), 9); // Assuming linvec outputs a 1x9 vector, prepare the output matrix
    
    for (int i = 0; i < projs1.cols(); ++i) {
        Eigen::VectorXd a = projs1.col(i); // Get the i-th column (equivalent to i-th transposed row)
        Eigen::VectorXd b = projs2.col(i); // Ditto for projs2
        coeffs.row(i) = linvec(a, b).transpose(); // Apply linvec and store the result in coeffs
    }

    return coeffs;
}
