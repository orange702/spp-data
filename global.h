#pragma once

#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

// WGS-84椭球模型及物理常量
const double PI = 3.141592653589793;
const double WGS84_A = 6378137.0;              // WGS-84 椭球长半轴
const double WGS84_F = 1.0 / 298.257223563;    // WGS-84 椭球扁率
const double WGS84_E2 = WGS84_F * (2.0 - WGS84_F); // WGS-84 第一偏心率的平方
const double C_LIGHT = 299792458.0;            // 真空中光速 (m/s)
const double GM_EARTH = 3.986005E14;           // 地球引力常数 (m^3/s^2)
const double OMEGA_EARTH = 7.2921151467E-5;   // 地球自转角速度 (rad/s)

// GPS星历数据结构
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

// 单个历元的观测数据
struct Observation {
    double time; // 观测时刻 (s) in GPS week
    double C1C;  // C1C 伪距观测值
    // 可在此处添加其他观测量 L1C, D1C等
};

// 测站信息
struct Station {
    std::string name;
    Eigen::Vector3d approx_pos; // 测站近似坐标 XYZ
    double cutoff_angle;       // 截止高度角 (degrees)
};

// 卫星计算结果
struct SatResult {
    Eigen::Vector3d pos; // 卫星位置 XYZ
    double clock_error;  // 卫星钟差 (s)
    double elevation;    // 高度角 (rad)
    double azimuth;      // 方位角 (rad)
};