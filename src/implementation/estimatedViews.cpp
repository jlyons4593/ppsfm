#include "estimatedViews.h"
#include "dataStructures.hpp"
#include "options.h"
#include <iostream>

EstimatedViews::EstimatedViews(DataStructures::SfMData data,Eigen::MatrixXd points, Eigen::VectorXi index_points,int new_view, int num_fixed){
    int num_points = index_points.size();

    Eigen::Vector2i new_view_idx(new_view, new_view+1);

    // Initialize LLS system
    sys = Eigen::MatrixXd::Zero(2 * num_points, 12);

    // Initialize linear equality constraint (fixes sum of the fixed projective depths)
    con = Eigen::RowVectorXd::Zero(12);

    // Initialize positivity constraint (projective depths)
    pos = Eigen::MatrixXd::Zero(num_points, 12);

    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(3, 12);


    for (int i = 0; i < num_points; ++i) {
    // for (int i = 0; i < 1; ++i) {

        Eigen::Vector3i point_idx(3*index_points(i), 3*index_points(i)+1, 3*index_points(i)+2);
        g.block<1, 4>(0, 0) = points.col(index_points[i]);
        g.block<1, 4>(1, 4) = points.col(index_points[i]);
        g.block<1, 4>(2, 8) = points.col(index_points[i]);
    
        Eigen::Matrix<double, 2, 3> selected_data;
        for (int k =0; k<new_view_idx.size(); k++){
            for (int j =0; j<point_idx.size(); j++){
                selected_data(k,j) = data.cost_function_data(new_view_idx(k), point_idx(j));
            }
        }
        sys.block<2, 12>(2 * i, 0) = selected_data* g;
        Eigen::Matrix<double, 1, 3> selected_pinv_meas;
        selected_pinv_meas.segment<3>(0) = data.pseudo_inverse_measurements.row(new_view).segment<3>(point_idx(0));
        pos.row(i) = selected_pinv_meas * g;
        if (i <= num_fixed) {
            con+= pos.row(i);
        }
        
    }
    con = con/num_fixed;

    Eigen::MatrixXd A = sys.rightCols(sys.cols() - 1) - sys.col(0) * con.segment(1, con.size() - 1) / con(0);
    Eigen::VectorXd b = -sys.col(0) / con(0);
    // Perform QR decomposition with column pivoting
    // Using ColPivHouseHolder which is more similar but doesnt provide the economy option that householder provides
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q_full = qr.matrixQ();
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(A.cols(), A.cols()).triangularView<Eigen::Upper>();
    Eigen::VectorXi p = qr.colsPermutation().indices();
    Eigen::MatrixXd Q = Q_full.leftCols(A.cols());

    // Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    // Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(A.rows(), A.cols());
    // Eigen::MatrixXd R = qr.matrixQR().topLeftCorner(A.cols(), A.cols()).triangularView<Eigen::Upper>();
    // Eigen::VectorXi p(A.cols());
    // p.setLinSpaced(A.cols(), 0, A.cols() - 1);

    bool hasSmallDiagonalElement = false;
    for (int i = 0; i < R.diagonal().size(); ++i) {
        if (std::abs(R.diagonal()(i)) < Options::RANK_TOLERANCE) {
            hasSmallDiagonalElement = true;
            break;
        }
    }
    if (hasSmallDiagonalElement){
        Eigen::VectorXd estimation; 
    }
    else{
        Eigen::MatrixXd d = Q.transpose()*b;
        Eigen::MatrixXd x = R.triangularView<Eigen::Upper>().solve(d);
        // std::cout<<p<<std::endl;
        // std::cout<<x<<std::endl;
        Eigen::VectorXi pp = Eigen::VectorXi::LinSpaced(11, 0, 11);
        // std::cout<<pp<<std::endl;
        // std::cout<<std::endl;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(p);
        pp = perm * pp;
        Eigen::VectorXd y = perm * x;

        estim.resize(x.size() + 1);
        double scalar = 1.0 / con(0);
        Eigen::VectorXd rowVec = con.segment(1, con.size() - 1);
        double dotProduct = rowVec.transpose().dot(y);
        double result = scalar -dotProduct/ con(0);
        estim(0) = result;
        estim.segment(1, y.size()) = y;
        // std::cout<<estim<<std::endl;

    }
}
