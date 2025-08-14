#pragma once

#include "rinex.h"

class SPPSolver {
public:
    // solve 函数的返回类型是 void
    void solve(const Station& station,
        const std::map<double, std::map<int, Observation>>& obs_data,
        const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
        const std::vector<double>& ion_params);

    // 【新增】calculateRMS 函数的声明
    void calculateRMS(const std::vector<Eigen::Vector3d>& positions, const Eigen::Vector3d& ref_pos);


private:
    // 为单个历元执行SPP解算
    bool processEpoch(double time,
        const std::map<int, Observation>& epoch_obs,
        const Station& station,
        const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
        const std::vector<double>& ion_params,
        Eigen::Vector4d& solution,
        int& valid_sats_count); // <--- 新增

    // 计算卫星位置和钟差
    bool calculateSatPosition(double obs_time, double pseudorange, const GPSEphemeris& eph, SatResult& result);

    // 计算卫星高度角和方位角
    void calculateAzimuthElevation(const Eigen::Vector3d& rec_pos, const Eigen::Vector3d& sat_pos, double& az, double& el);

    // 从星历库中选择最合适的星历
    const GPSEphemeris* selectEphemeris(double time, int prn, const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris);
};