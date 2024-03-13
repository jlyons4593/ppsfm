#include <Eigen/Dense>
#pragma once
class PyramidalVisibilityScore {
private:
    int level;
    int width;
    int height;
    Eigen::MatrixXi proj_count;
    Eigen::RowVectorXi width_range;
    Eigen::RowVectorXd height_range;
    Eigen::RowVectorXd dimension_range;

    std::pair<Eigen::RowVectorXi,Eigen::RowVectorXi> visibleCell(Eigen::MatrixXd projections);
    
    void addProjections(const Eigen::MatrixXd projections); 
    void removeProjections(const Eigen::MatrixXi& projs); 

    Eigen::MatrixXi downsampleAndCombine(const Eigen::MatrixXi& visibility); 

    double maxScore(); 

public:
    PyramidalVisibilityScore(int image_width, int image_height, int score_level = -1, const Eigen::MatrixXd& projections = Eigen::MatrixXd());

    double computeScore(bool normalized = false);
};

