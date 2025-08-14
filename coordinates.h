#pragma once
#include "global.h"

// XYZ to BLH (Latitude, Longitude, Height)
Eigen::Vector3d xyz2blh(const Eigen::Vector3d& xyz);

// XYZ to NEU (North, East, Up)
Eigen::Vector3d xyz2neu(const Eigen::Vector3d& xyz, const Eigen::Vector3d& ref_xyz);