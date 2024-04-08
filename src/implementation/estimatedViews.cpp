#include "estimatedViews.h"
#include "dataStructures.hpp"
#include "options.h"
#include <exception>
#include <iostream>

EstimatedViews::EstimatedViews(DataStructures::SfMData data,Eigen::MatrixXd points, Eigen::VectorXi index_points,int new_view, int num_fixed){
    int num_points = index_points.size();

    // std::cout<<new_view<<std::endl;
    Eigen::Vector2i new_view_idx(2*new_view, (2*new_view)+1);
    // std::cout<<new_view_idx<<std::endl;

    // Initialize LLS system
    sys = Eigen::MatrixXd::Zero(2 * num_points, 12);

    // Initialize linear equality constraint (fixes sum of the fixed projective depths)
    con = Eigen::RowVectorXd::Zero(12);

    // Initialize positivity constraint (projective depths)
    pos = Eigen::MatrixXd::Zero(num_points, 12);

    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(3, 12);
    


    // auto start = std::chrono::high_resolution_clock::now();
    Eigen::Matrix<double, 1, 3> selected_pinv_meas;
    Eigen::Vector3i point_idx(3);
    Eigen::Matrix<double, 2, 3> selected_data;
    for (int i = 0; i < num_points; ++i) {
        point_idx(0) = 3*index_points(i);
        point_idx(1) = 3*index_points(i)+1;
        point_idx(2) = 3*index_points(i)+2;
        g.block<1, 4>(0, 0) = points.col(index_points[i]);
        g.block<1, 4>(1, 4) = points.col(index_points[i]);
        g.block<1, 4>(2, 8) = points.col(index_points[i]);
    
        // look for a way to speed this up
        selected_data(0,0) = data.cost_function_data(new_view_idx(0), point_idx(0));
        selected_data(1,0) = data.cost_function_data(new_view_idx(1), point_idx(0));
        selected_data(0,1) = data.cost_function_data(new_view_idx(0), point_idx(1));
        selected_data(1,1) = data.cost_function_data(new_view_idx(1), point_idx(1));
        selected_data(0,2) = data.cost_function_data(new_view_idx(0), point_idx(2));
        selected_data(1,2) = data.cost_function_data(new_view_idx(1), point_idx(2));
        sys.block<2, 12>(2 * i, 0) = selected_data* g;

        selected_pinv_meas.block(0,0,1,3) = data.pseudo_inverse_measurements.block(new_view, point_idx(0), 1,3);
        pos.row(i) = selected_pinv_meas * g;
        if (i < num_fixed) {
            con+=pos.row(i);
        }
    }

    // std::cout<<sys<<std::endl;
    con = con/num_fixed;

    Eigen::MatrixXd A = sys.rightCols(sys.cols() - 1) - sys.col(0) * con.segment(1, con.size() - 1) / con(0);
    Eigen::VectorXd b = -sys.col(0) / con(0);
    // std::cout<<"B"<<std::endl;
    // std::cout<<b<<std::endl;
    // Perform QR decomposition with column pivoting
    // Using ColPivHouseHolder which is more similar but doesnt provide the economy option that householder provides

    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    // Eigen::MatrixXd QR = qr.matrixQR();
    // auto end = std::chrono::high_resolution_clock::now();
    // bool timer = true;
    // if (timer) {
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //     std::cout << "Loop execution in " << duration.count() << " milliseconds" << std::endl;
    // }

    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    // THis line is till very slow
    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    // Eigen::MatrixXd Q = qr.matrixQ();
    // Eigen::MatrixXd R = qr.matrixR().topLeftCorner(A.cols(), A.cols()).triangularView<Eigen::Upper>();
    // Eigen::VectorXi p = qr.colsPermutation().indices();
    // Eigen::MatrixXd Q = Q_full.leftCols(A.cols());
    // DataStructures::printColsRows(Q, "Q");

    // Faster version using householderQ
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd QR = qr.matrixQR();
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(A.cols(), A.cols()).triangularView<Eigen::Upper>();
    Eigen::VectorXi p = qr.colsPermutation().indices();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    if(duration.count()>1){
        std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    
    }
    // Eigen::MatrixXd Q = Q_full.leftCols(A.cols());
    // DataStructures::printColsRows(Q, "Q");
    Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(A.rows(),A.cols()));
    thinQ = qr.householderQ()*thinQ;
        // Eigen::MatrixXd d = thinQ.transpose()*b;


    bool hasSmallDiagonalElement = false;
    for (int i = 0; i < R.diagonal().size(); ++i) {
        if (std::abs(R.diagonal()(i)) < 1e-6) {
            hasSmallDiagonalElement = true;
            break;
        }
    }
    if (hasSmallDiagonalElement){
        Eigen::VectorXd estimation; 
    }
    else{
        Eigen::MatrixXd d = thinQ.transpose()*b;
        // Eigen::MatrixXd d = Q.transpose()*b;
    // auto end = std::chrono::high_resolution_clock::now();
    // bool timer = true;
    // if (timer) {
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //     std::cout << "Loop execution in " << duration.count() << " milliseconds" << std::endl;
    // }
        Eigen::MatrixXd x = R.triangularView<Eigen::Upper>().solve(d);
        // std::cout<<p<<std::endl;
        // std::cout<<x<<std::endl;
        Eigen::VectorXi pp = Eigen::VectorXi::LinSpaced(11, 0, 11);
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

    }
}
