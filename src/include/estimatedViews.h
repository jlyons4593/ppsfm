
#include "dataStructures.hpp"

class EstimatedViews{
  private:

  public:

    Eigen::MatrixXd sys;  

    // Initialize linear equality constraint (fixes sum of the fixed projective depths)
    Eigen::RowVectorXd con; 

    // Initialize positivity constraint (projective depths)
    Eigen::MatrixXd pos; 

    Eigen::VectorXd estim;
    
    EstimatedViews(DataStructures::SfMData data,Eigen::MatrixXd points, Eigen::VectorXi index_points,int new_view, int num_fixed);

};
