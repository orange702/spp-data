#pragma once
#include "global.h"

// Klobuchar 电离层模型
double klobucharModel(const Eigen::Vector3d& rec_pos_blh, double az, double el, double gpst, const std::vector<double>& ion_params);

// Saastamoinen 对流层模型
double saastamoinenModel(double el, double rec_height);