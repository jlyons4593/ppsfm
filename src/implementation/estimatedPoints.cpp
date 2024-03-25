
#include "estimatedPoints.h"
#include "options.h"
    
EstimatedPoints::EstimatedPoints(DataStructures::SfMData& data,Eigen::MatrixXd& cameras, Eigen::VectorXi& idx_views,int new_point, int num_fixed){

    int num_views = idx_views.size();

    Eigen::Vector3i new_point_idx(3*new_point, (3*new_point)+1, (3*new_point)+2);

    //  Initialize LLS system
    sys = Eigen::MatrixXd::Zero(2 * num_views, 4);

    // Initialize linear equality constraint (fixes sum of the fixed projective depths)
    con = Eigen::RowVectorXd::Zero(4);

    // Initialize positivity constraint (projective depths)
    pos = Eigen::MatrixXd::Zero(num_views, 4);

    
    for (int j = 0; j < num_views; ++j) {
    // for (int j = 3; j < 4; ++j) {
        Eigen::Vector3i view_idx(3*idx_views(j), (3*idx_views(j))+1, (3*idx_views(j))+2);
        Eigen::MatrixXd selected_data = data.cost_function_data.block(2 * idx_views(j), new_point_idx(0), 2, 3);
        // std::cout<<selected_data<<std::endl;
        Eigen::MatrixXd selected_cameras = cameras.block(view_idx(0), 0, 3, cameras.cols());
        // std::cout<<selected_cameras<<std::endl;

        sys.block(j, 0, 2,sys.cols())= selected_data*selected_cameras;


        Eigen::MatrixXd selected_pinv_meas = data.pseudo_inverse_measurements.block(idx_views(j),0,1,3);
        // std::cout<<selected_pinv_meas<<std::endl;
        
        pos.row(j) = selected_pinv_meas*selected_cameras;

        if(j<num_fixed) {
            con.row(0) = con.row(0) + pos.row(j);
        }
    }
    

    con = con/num_fixed;

    Eigen::MatrixXd A = sys.rightCols(sys.cols() - 1) - sys.col(0) * con.segment(1, con.size() - 1) / con(0);
    Eigen::VectorXd b = -sys.col(0) / con(0);

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q_full = qr.matrixQ();
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(A.cols(), A.cols()).triangularView<Eigen::Upper>();
    Eigen::VectorXi p = qr.colsPermutation().indices();
    Eigen::MatrixXd Q = Q_full.leftCols(A.cols());
            

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
        Eigen::VectorXi pp = Eigen::VectorXi::LinSpaced(3, 0, 3);
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

