// coordinates.cpp (The Absolutely, Finally, Correct Version)

#include "coordinates.h"
#include <cmath>

// XYZ to BLH (Latitude, Longitude, Height)
// 这个函数的实现是正确的，无需修改
Eigen::Vector3d xyz2blh(const Eigen::Vector3d& xyz) {
    double x = xyz[0], y = xyz[1], z = xyz[2];
    double p = sqrt(x * x + y * y);
    double lat, lon, height;

    // 迭代法求解纬度，精度很高
    lat = atan2(z, p * (1.0 - WGS84_E2));
    double N;
    for (int i = 0; i < 5; ++i) {
        N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin(lat) * sin(lat));
        height = p / cos(lat) - N;
        lat = atan2(z, p * (1.0 - WGS84_E2 * N / (N + height)));
    }
    lon = atan2(y, x);

    // 重新计算最后一次的卯酉圈曲率半径和高程
    N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin(lat) * sin(lat));
    height = p / cos(lat) - N;

    return { lat, lon, height };
}


// XYZ to NEU (North, East, Up)
// 【最终修正版】修正了旋转矩阵的错误
Eigen::Vector3d xyz2neu(const Eigen::Vector3d& xyz, const Eigen::Vector3d& ref_xyz) {
    // 1. 计算参考点的经纬度
    Eigen::Vector3d ref_blh = xyz2blh(ref_xyz);
    double lat = ref_blh[0];
    double lon = ref_blh[1];

    // 2. 计算待求点相对于参考点的偏差向量
    Eigen::Vector3d d_xyz = xyz - ref_xyz;

    // 3. 构建从XYZ到NEU的正确旋转矩阵 R
    Eigen::Matrix3d R;
    R(0, 0) = -sin(lat) * cos(lon); R(0, 1) = -sin(lat) * sin(lon); R(0, 2) = cos(lat);        // North row
    R(1, 0) = -sin(lon);            R(1, 1) = cos(lon);             R(1, 2) = 0;                  // East row
    R(2, 0) = cos(lat) * cos(lon);  R(2, 1) = cos(lat) * sin(lon);  R(2, 2) = sin(lat);        // Up row

    // 4. 执行旋转
    return R * d_xyz;
}