// coordinates.cpp (The Absolutely, Finally, Correct Version)

#include "coordinates.h"
#include <cmath>

// XYZ to BLH (Latitude, Longitude, Height)
// ���������ʵ������ȷ�ģ������޸�
Eigen::Vector3d xyz2blh(const Eigen::Vector3d& xyz) {
    double x = xyz[0], y = xyz[1], z = xyz[2];
    double p = sqrt(x * x + y * y);
    double lat, lon, height;

    // ���������γ�ȣ����Ⱥܸ�
    lat = atan2(z, p * (1.0 - WGS84_E2));
    double N;
    for (int i = 0; i < 5; ++i) {
        N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin(lat) * sin(lat));
        height = p / cos(lat) - N;
        lat = atan2(z, p * (1.0 - WGS84_E2 * N / (N + height)));
    }
    lon = atan2(y, x);

    // ���¼������һ�ε�î��Ȧ���ʰ뾶�͸߳�
    N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin(lat) * sin(lat));
    height = p / cos(lat) - N;

    return { lat, lon, height };
}


// XYZ to NEU (North, East, Up)
// �����������桿��������ת����Ĵ���
Eigen::Vector3d xyz2neu(const Eigen::Vector3d& xyz, const Eigen::Vector3d& ref_xyz) {
    // 1. ����ο���ľ�γ��
    Eigen::Vector3d ref_blh = xyz2blh(ref_xyz);
    double lat = ref_blh[0];
    double lon = ref_blh[1];

    // 2. ������������ڲο����ƫ������
    Eigen::Vector3d d_xyz = xyz - ref_xyz;

    // 3. ������XYZ��NEU����ȷ��ת���� R
    Eigen::Matrix3d R;
    R(0, 0) = -sin(lat) * cos(lon); R(0, 1) = -sin(lat) * sin(lon); R(0, 2) = cos(lat);        // North row
    R(1, 0) = -sin(lon);            R(1, 1) = cos(lon);             R(1, 2) = 0;                  // East row
    R(2, 0) = cos(lat) * cos(lon);  R(2, 1) = cos(lat) * sin(lon);  R(2, 2) = sin(lat);        // Up row

    // 4. ִ����ת
    return R * d_xyz;
}