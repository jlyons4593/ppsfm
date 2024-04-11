#include "Helper.hpp"
#include <iostream>
#include <random>
#include <unordered_set>

namespace Helper{
    Eigen::VectorXi random_subset(const Eigen::VectorXi& complete_set, int num_rejected, int num_sample) {
        int n = complete_set.size();
        static std::mt19937 gen(std::random_device{}()); // Random number generator

        if (num_rejected > 0) {

            std::uniform_int_distribution<int> dist(0, complete_set.size() - num_rejected);

            // Generate a random integer
            int random_value = dist(gen);

            // Calculate non_rejected
            int non_rejected = num_rejected + random_value;
            // Eigen::VectorXi temp_set(num_sample - 1);
            // std::uniform_int_distribution<> dis_sample(0, n - 3);
            // for (int i = 0; i < num_sample - 1; ++i) {
            //     int idx = dis_sample(gen);
            //     if (idx >= non_rejected)
            //         ++idx;
            //     temp_set(i) = idx;
            // }
            Eigen::VectorXi temp_set(num_sample - 1);
            std::unordered_set<int> selected_indices;
            std::uniform_int_distribution<> dis_sample(0, n - 3);

            int count = 0;
            while (count < num_sample - 1) {
                int idx = dis_sample(gen);
                if (idx >= non_rejected)
                    ++idx;

                if (selected_indices.find(idx) == selected_indices.end()) {
                    selected_indices.insert(idx);
                    temp_set(count) = idx;
                    ++count;
                }
            }
            // temp_set(temp_set >= non_rejected) = temp_set(temp_set >= non_rejected) + 1;
            // Not needed since we already handled this in the loop above

            // test_set = horzcat(temp_set, non_rejected);
            Eigen::VectorXi test_set(num_sample);
            test_set.head(num_sample - 1) = temp_set;
            test_set(num_sample - 1) = non_rejected-1;

            return test_set;
        } else {
            Eigen::VectorXi indices = Eigen::VectorXi::LinSpaced(n, 0, n - 1);

            std::shuffle(indices.data(), indices.data() + n, gen);

            // Take the first num_sample elements as the test set
            return indices.head(num_sample);
        }
    }

} 
