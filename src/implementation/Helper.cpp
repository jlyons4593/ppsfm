#include "Helper.hpp"
#include <iostream>
#include <random>

namespace Helper{
    Eigen::VectorXi random_subset(const Eigen::VectorXi& complete_set, int num_rejected, int num_sample) {
        int n = complete_set.size();
        static std::mt19937 gen(std::random_device{}()); // Random number generator

        if (num_rejected > 0) {
            std::uniform_int_distribution<> dis(num_rejected, n - 1);
            int non_rejected = dis(gen);

            // temp_set = randperm(length(complete_set)-1, num_sample-1);
            Eigen::VectorXi temp_set(num_sample - 1);
            std::uniform_int_distribution<> dis_sample(0, n - 2);
            for (int i = 0; i < num_sample - 1; ++i) {
                int idx = dis_sample(gen);
                if (idx >= non_rejected)
                    ++idx;
                temp_set(i) = idx;
            }

            // temp_set(temp_set >= non_rejected) = temp_set(temp_set >= non_rejected) + 1;
            // Not needed since we already handled this in the loop above

            // test_set = horzcat(temp_set, non_rejected);
            Eigen::VectorXi test_set(num_sample);
            test_set.head(num_sample - 1) = temp_set;
            test_set(num_sample - 1) = non_rejected;

            return test_set;
        } else {
            int n = complete_set.size();
            std::random_device rd;
            std::mt19937 gen(rd());

            Eigen::VectorXi indices = Eigen::VectorXi::LinSpaced(n, 0, n - 1);
            std::shuffle(indices.data(), indices.data() + n, gen);

            Eigen::VectorXi test_set(num_sample);
            test_set.head(num_sample) = indices.head(num_sample);

            return test_set;
        }
    }

} 
