#pragma once
#include "global.h"

// Klobuchar �����ģ��
double klobucharModel(const Eigen::Vector3d& rec_pos_blh, double az, double el, double gpst, const std::vector<double>& ion_params);

// Saastamoinen ������ģ��
double saastamoinenModel(double el, double rec_height);