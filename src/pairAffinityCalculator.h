#include "stage.h"
class PairAffinityStage : public Stage {
public:
    void process(std::shared_ptr<InputData> input) override {
        auto affinityInput = std::dynamic_pointer_cast<AffinityInputData>(input);
        if (affinityInput) {
            // Process affinityInput->pairs
        }
    }
};
