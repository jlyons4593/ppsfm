
#include "dataStructures.hpp"
#include <iostream>
class Diagnosis{

  public:
    Diagnosis(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& inliers, DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables, int iteration);

};
