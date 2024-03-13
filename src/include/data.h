
#include <Eigen/Sparse>
class InputData {
public:
    virtual ~InputData() = default; // Ensure polymorphic deletion
};

class CleaningInputData : public InputData {
public:
    Eigen::SparseMatrix<double> measurements;
    // Additional specific attributes for data cleaning
};

class AffinityInputData : public InputData {
public:
    // Additional specific attributes for affinity generation
};

