#pragma once

#include "rinex.h"

class SPPSolver {
public:
    // solve �����ķ��������� void
    void solve(const Station& station,
        const std::map<double, std::map<int, Observation>>& obs_data,
        const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
        const std::vector<double>& ion_params);

    // ��������calculateRMS ����������
    void calculateRMS(const std::vector<Eigen::Vector3d>& positions, const Eigen::Vector3d& ref_pos);


private:
    // Ϊ������Ԫִ��SPP����
    bool processEpoch(double time,
        const std::map<int, Observation>& epoch_obs,
        const Station& station,
        const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
        const std::vector<double>& ion_params,
        Eigen::Vector4d& solution,
        int& valid_sats_count); // <--- ����

    // ��������λ�ú��Ӳ�
    bool calculateSatPosition(double obs_time, double pseudorange, const GPSEphemeris& eph, SatResult& result);

    // �������Ǹ߶ȽǺͷ�λ��
    void calculateAzimuthElevation(const Eigen::Vector3d& rec_pos, const Eigen::Vector3d& sat_pos, double& az, double& el);

    // ����������ѡ������ʵ�����
    const GPSEphemeris* selectEphemeris(double time, int prn, const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris);
};