#include "dataStructures.hpp"
#pragma once
class Refinement{
  private:
    DataStructures::SfMData data;
    DataStructures::ComputedCameraPoints camera_variables;
  public:
    Refinement(DataStructures::SfMData& data, DataStructures::ComputedCameraPoints& camera_variables , int init_refine, int last_path, bool start_cameras, int type);


};
