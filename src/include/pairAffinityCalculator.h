#include "dataStructures.hpp"
#include "pyramidalVisibilityScore.h"
#include <Eigen/Dense>
class PairAffinityCalculator{
private:
  DataStructures::ViewpairAffinity pairAffinity;

public:
  PairAffinityCalculator() {}
  DataStructures::ViewpairAffinity getPairAffinity();
  void process(Eigen::MatrixXd image_measurements, Eigen::MatrixXd visible,
               Eigen::MatrixXd image_sizes);
   
};
