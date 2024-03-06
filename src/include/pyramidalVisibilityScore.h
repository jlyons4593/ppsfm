#include <Eigen/Dense>
#include "dataStructures.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

class PyramidalVisibilityScore {
private:
    int level;
    int width;
    int height;
    Eigen::MatrixXi proj_count;
    Eigen::RowVectorXi width_range;
    Eigen::RowVectorXd height_range;
    Eigen::RowVectorXd dimension_range;
    std::pair<Eigen::RowVectorXi,Eigen::RowVectorXi> visibleCell(Eigen::MatrixXd projections){
        // 
        Eigen::RowVectorXi idx_width = Eigen::RowVectorXi::Ones(projections.cols());
        Eigen::RowVectorXi idx_height = Eigen::RowVectorXi::Ones(projections.cols());

        // // std::cout<<idx_width.size()<<" "<<idx_height.size()<<std::endl;
        // // for (int k = 0; k < dimension_range.size(); ++k) {
        //
        for (int k = 0; k < dimension_range.size(); ++k) {
            Eigen::RowVectorXi idx_middle_width = idx_width + Eigen::RowVectorXi::Constant(projections.cols(), dimension_range(k));
            Eigen::RowVectorXi idx_middle_height = idx_height + Eigen::RowVectorXi::Constant(projections.cols(), dimension_range(k));

            //Width section
            Eigen::VectorXd comparisonVector(idx_middle_width.size());
            for (int i = 0; i < idx_middle_width.size(); ++i) {
                // Adjust for 0-based indexing by subtracting 1 from idx_middle(i)
                comparisonVector(i) = width_range(idx_middle_width(i) );
            }

            Eigen::Array<bool, Eigen::Dynamic, 1> idx_right(projections.cols());
            for (int i = 0; i < projections.cols(); ++i) {

                idx_right(i) = projections(0, i) > comparisonVector(i);
            }
             for (int i = 0; i < idx_right.size(); ++i) {
                if (idx_right(i)) { // If the condition is true
                    idx_width(i) = idx_middle_width(i);
                }
            }
            //Height section

            Eigen::VectorXd comparisonVector2(idx_middle_height.size());
            for (int i = 0; i < idx_middle_height.size(); ++i) {
                // Adjust for 0-based indexing by subtracting 1 from idx_middle(i)
                comparisonVector2(i) = height_range(idx_middle_height(i));
            }

            Eigen::Array<bool, Eigen::Dynamic, 1> idx_height_right(projections.cols());
            for (int i = 0; i < projections.cols(); ++i) {
                
                idx_height_right(i) = projections(1, i) > comparisonVector2(i);
            }
             for (int i = 0; i < idx_height_right.size(); ++i) {
                if (idx_height_right(i)) { // If the condition is true
                    idx_height(i) = idx_middle_height(i);
                }
            }
        }

      return {idx_width,idx_height};
    }
public:
    PyramidalVisibilityScore(int image_width, int image_height, int score_level = -1, const Eigen::MatrixXd& projections = Eigen::MatrixXd())
        : width(image_width), height(image_height), level(score_level) {

        int dimension = std::pow(2, level); 

        proj_count = Eigen::MatrixXi::Zero(dimension,dimension);
        width_range = Eigen::VectorXi::LinSpaced(std::pow(2,level), 0, width);
        height_range= Eigen::VectorXd::LinSpaced(std::pow(2,level)+1, 0, height);
        

        dimension_range.resize(level);
        for (int i = 0; i<level; ++i) {
            // index to -1 to make sure it is in bounds
            // starting from level so that it reads big to small
            dimension_range(level-i-1) = std::pow(2, i);
        }

        addProjections(projections);

    }


    void addProjections(const Eigen::MatrixXd projections) {
        if(projections.cols()>2){
            auto [idx_width, idx_height] = visibleCell(projections);
            // std::cout<<"idx width: "<<idx_width<<std::endl;
            // std::cout<<"idx width size: "<<idx_width.size()<<std::endl;
            // std::cout<<"idx height: "<<idx_height<<std::endl;
            // std::cout<<"idx height size: "<<idx_height.size()<<std::endl;
           for(int i=0; i<idx_width.size();i++) {

                proj_count(idx_height(i), idx_width(i))++;
            }
        }
        // std::cout<<proj_count<<std::endl;
    }


    void removeProjections(const Eigen::MatrixXi& projs) {
        // Similar to addProjections, but decrementing the proj_count for each projection
    }
    Eigen::MatrixXi downsampleAndCombine(const Eigen::MatrixXi& visibility) {
        // Assuming visibility is a binary matrix (elements are 0 or 1)
        // Create a new matrix half the size of 'visibility'
        Eigen::MatrixXi reduced(visibility.rows() / 2, visibility.cols() / 2);

        for (int i = 0; i < reduced.rows(); ++i) {
            for (int j = 0; j < reduced.cols(); ++j) {
                // Perform logical OR on four adjacent elements from 'visibility'
                reduced(i, j) = visibility(2 * i, 2 * j) || visibility(2 * i + 1, 2 * j) ||
                    visibility(2 * i, 2 * j + 1) || visibility(2 * i + 1, 2 * j + 1);
            }
        }

        return reduced;
    }
    double computeScore(bool normalized = false) {
        // Implement the score computation as in MATLAB, considering the Eigen way of handling matrices
        // and the logical operations for visibility
        int score=0;
        // there should be a check to see if the projection_count is empty
        if(level>0){
            Eigen::MatrixXi visibility = (proj_count.array() > 0).cast<int>();
            // std::cout<<visibility<<std::endl;
            for(int k = 0; k<dimension_range.size(); k++){
                int nonZeroCount=0;
                for(int i=0; i< visibility.cols(); i++){
                    for(int j=0; j< visibility.rows(); j++){
                        if (visibility(i,j)>0){
                            nonZeroCount++;
                        }
                    }
                }
                score = score + nonZeroCount * dimension_range(k);
                visibility = downsampleAndCombine(visibility);
            }
            //visib filled

        }
        return score;
    }

    double maxScore() {
        int score = 0;
        if(level>0){
            // Perform the operations
            Eigen::VectorXd result = dimension_range.array() * 2; // Multiply each element by 2
            result = result.array().square(); // Square each element
            double score = (dimension_range.array() * result.array()).sum(); // Multiply element-wise with original and sum
        }
        return score;
        // Implement the maximum score computation
    }


    // Private methods if needed
};

