#include "dataStructures.hpp"
#include "pyramidalVisibilityScore.h"
#pragma once
class FactorCompletion{
private:
    DataStructures::SfMData data;
    DataStructures::ComputedCameraPoints camera_variables;
    DataStructures::ViewpairAffinity pair_affinity;
    int last_path;

    Eigen::VectorXi rejected_views;
    Eigen::VectorXi rejected_points;
    Eigen::MatrixXd centers;
    Eigen::MatrixXd image_sizes;
    std::vector<PyramidalVisibilityScore*> pvs_scores;
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> inliers;
public:

    FactorCompletion(DataStructures::SfMData data,   
    DataStructures::ComputedCameraPoints camera_variables,
    DataStructures::ViewpairAffinity pair_affinity, Eigen::MatrixXd image_sizes, Eigen::MatrixXd centers): data(data), camera_variables(camera_variables), pair_affinity(pair_affinity), image_sizes(image_sizes), centers(centers) {}

    void process();

    void check_expand_init(int num_points, int num_views);

    void init_pvs(int num_views);
    
    std::pair<Eigen::RowVectorXd, std::vector<int>> search_eligible_views(Eigen::VectorXi thresholds, Eigen::VectorXi rejected_views);

    Eigen::VectorXd search_eligible_points(int required_visible_eligibility, Eigen::VectorXi rejected_points);

    int try_adding_views(std::vector<int> eligibles, int level_views);
};


