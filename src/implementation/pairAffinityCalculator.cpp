#include "pairAffinityCalculator.h"
#include "logger.h"
#include "options.h"
#include <numeric>

  void PairAffinityCalculator::process(Eigen::MatrixXd image_measurements, Eigen::MatrixXd visible,
               Eigen::MatrixXd image_sizes) {

    // Logger::logSection("Pair Affinity");

    int num_views = visible.rows();
    // std::cout<<image_sizes<<std::endl;

    // Convert 'visible' to a logical matrix
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> logicalVisible = (visible.array() > 0).cast<bool>();
    int num_pair = (num_views * (num_views - 1)) / 2;
    Eigen::MatrixXd view_pairs = Eigen::MatrixXd::Zero(num_pair, 2);
    // std::cout<<view_pairs.size()<<std::endl;
    Eigen::VectorXd affinity = Eigen::VectorXd::Zero(num_pair);

    int offset = 0;

    // Logger::logSubsection("Processing Views Using Pyramidal Visibility Score");

    for (int first_view = 0; first_view<num_views -1; ++first_view) {
// #pragma omp parallel for
      for (int second_view = first_view + 1; second_view < num_views;
      ++second_view) {


        // Logical AND to find common points
        Eigen::ArrayXi common_points = (logicalVisible.row(first_view).array() && logicalVisible.row(second_view).array()).cast<int>();
        // std::cout<<common_points<<std::endl;

        int num_common = common_points.count(); // Count non-zeros

        if (num_common>Options::MIN_COMMON_INIT){
          view_pairs(offset, 0) = first_view;
          view_pairs(offset, 1) = second_view;


          Eigen::MatrixXd image_measurements_chunk1(2,num_common);
          Eigen::MatrixXd image_measurements_chunk2(2,num_common);

          int count=0;
          for(int i =0; i<common_points.size(); i++){

            if(common_points(i)>0){

              image_measurements_chunk1(0, count) = image_measurements(first_view*3, i); 
              image_measurements_chunk1(1, count) = image_measurements(first_view*3+1, i);
              image_measurements_chunk2(0, count) = image_measurements(second_view*3, i); 
              image_measurements_chunk2(1, count) = image_measurements(second_view*3+1, i);
              count++;
            }
          }

          PyramidalVisibilityScore pvs_first = PyramidalVisibilityScore(image_sizes(0), image_sizes(1), Options::SCORE_LEVEL, image_measurements_chunk1);
          PyramidalVisibilityScore pvs_second = PyramidalVisibilityScore(image_sizes(0), image_sizes(1), Options::SCORE_LEVEL, image_measurements_chunk2);


          int one = pvs_first.computeScore(false);
          int two = pvs_second.computeScore(false);
          // std::cout<<one<<":"<<two<<std::endl;
          affinity(offset) = one +two;  
          // std::cout<<affinity(offset)<<std::endl;

          offset++;

        }
      }
    }

    // Logger::logSubsection("Sorting View Pairs and Affinity");


    view_pairs = view_pairs.block(0, 0, offset, view_pairs.cols());
    affinity = affinity.block(0,0,offset,affinity.cols());
    // std::cout<<affinity.size()<<std::endl;
    // std::cout<<"view pair: "<<view_pairs<<std::endl;

    std::vector<int> indices(affinity.size());
    std::iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., affinity.size() - 1

    // Sort the indices based on comparing values in affinity
    std::sort(indices.begin(), indices.end(),
              [&affinity](int i1, int i2) { return affinity(i1) > affinity(i2); });

    // Create a sorted version of affinity using the sorted indices
    Eigen::VectorXd sortedAffinity(affinity.size());
    for (int i = 0; i < indices.size(); ++i) {
      sortedAffinity(i) = affinity(indices[i]);
    }

    Eigen::MatrixXd temp_view_pairs = Eigen::MatrixXd::Zero(view_pairs.rows(), view_pairs.cols());

    for (int i = 0; i < indices.size(); ++i) {
      temp_view_pairs.row(i) = view_pairs.row(indices[i]);
    }

    // Copy the sorted rows back to view_pairs
    view_pairs = temp_view_pairs;
    this->pairAffinity.view_pairs = view_pairs;
    this->pairAffinity.Affinity= sortedAffinity;
    // Logger::logSubsection("View Pairs and Affinity Stored as Member");

  } 

  DataStructures::ViewpairAffinity PairAffinityCalculator::getPairAffinity() {
    return this->pairAffinity;
  }
