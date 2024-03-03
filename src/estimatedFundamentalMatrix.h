
#include <Eigen/Dense>
#include <iostream>

class EstimatedFundamentalMatrix{
private:
    Eigen::MatrixXd fundamental_matrix;
    Eigen::VectorXi inliers;
    // std::pair<Eigen::MatrixXd,Eigen::VectorXi> ransac_estimate(){
    //
    // }

Eigen::VectorXd linvec(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
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
public:
    Eigen::MatrixXd getFundamentalMatrix(){

        return fundamental_matrix;
    }
    Eigen::VectorXi getInliers(){

        return inliers;
    }
    EstimatedFundamentalMatrix( const Eigen::MatrixXd& projs1, // Projections in the first image (2xN)
                               const Eigen::MatrixXd& projs2, // Corresponding projections in the second image (2xN)
                               double confidence,             // Confidence to stop robust estimation early in RANSAC
                               int max_iter,                  // Maximum number of iterations in RANSAC
                               double dist_thresh             // Distance threshold for computing inliers in RANSAC
                               )
    {

        
        // Eigen::MatrixXd projs1T = projs1.transpose();
        // Eigen::MatrixXd projs2T = projs2.transpose();
        // Eigen::MatrixXd combinedProjections(projs1T.rows(), projs1T.cols() + projs2T.cols());
        //
        // combinedProjections << projs1T, projs2T;
        // Currently only works for Ransac Method
        // std::pair<Eigen::MatrixXd, Eigen::VectorXi> fund_mat_inliers = ransac_estimate();

        
    }
};
