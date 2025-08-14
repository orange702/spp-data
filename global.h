#pragma once

#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

// WGS-84����ģ�ͼ�������
const double PI = 3.141592653589793;
const double WGS84_A = 6378137.0;              // WGS-84 ���򳤰���
const double WGS84_F = 1.0 / 298.257223563;    // WGS-84 �������
const double WGS84_E2 = WGS84_F * (2.0 - WGS84_F); // WGS-84 ��һƫ���ʵ�ƽ��
const double C_LIGHT = 299792458.0;            // ����й��� (m/s)
const double GM_EARTH = 3.986005E14;           // ������������ (m^3/s^2)
const double OMEGA_EARTH = 7.2921151467E-5;   // ������ת���ٶ� (rad/s)

// GPS�������ݽṹ
struct GPSEphemeris {
    int prn;
    double toc; // Clock data reference time (s) in GPS week
    double af0, af1, af2; // Clock correction parameters
    double iode, crs, delta_n, m0;
    double cuc, e, cus, sqrtA;
    double toe, cic, omega0, cis;
    double i0, crc, w, omega_dot;
    double idot;
    double tgd; // Group delay
    int week;   // GPS Week Number
    bool health;
};

// ������Ԫ�Ĺ۲�����
struct Observation {
    double time; // �۲�ʱ�� (s) in GPS week
    double C1C;  // C1C α��۲�ֵ
    // ���ڴ˴���������۲��� L1C, D1C��
};

// ��վ��Ϣ
struct Station {
    std::string name;
    Eigen::Vector3d approx_pos; // ��վ�������� XYZ
    double cutoff_angle;       // ��ֹ�߶Ƚ� (degrees)
};

// ���Ǽ�����
struct SatResult {
    Eigen::Vector3d pos; // ����λ�� XYZ
    double clock_error;  // �����Ӳ� (s)
    double elevation;    // �߶Ƚ� (rad)
    double azimuth;      // ��λ�� (rad)
};