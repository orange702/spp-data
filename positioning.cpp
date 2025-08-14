#include "positioning.h"
#include "models.h"
#include "coordinates.h"
#include <iostream>
#include <iomanip>

// ... solve() 和 calculateRMS() 函数完全不变 ...

void SPPSolver::solve(const Station& station,
    const std::map<double, std::map<int, Observation>>& obs_data,
    const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
    const std::vector<double>& ion_params) {

    std::vector<Eigen::Vector3d> solved_positions;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Epoch Time (s) |      X (m)      |      Y (m)      |      Z (m)      | clk (m) |  N (m)  |  E (m)  |  U (m)  | #Sats" << std::endl;
    std::cout << "----------------------------------------------------------------------------------------------------------------------" << std::endl;

    for (const auto& epoch_pair : obs_data) {
        double time = epoch_pair.first;
        const auto& epoch_obs = epoch_pair.second;

        Eigen::Vector4d solution;
        int valid_sats_count = 0;

        if (processEpoch(time, epoch_obs, station, gps_ephemeris, ion_params, solution, valid_sats_count)) {
            Eigen::Vector3d neu = xyz2neu({ solution[0], solution[1], solution[2] }, station.approx_pos);

            std::cout << std::setw(16) << time
                << std::setw(17) << solution[0]
                << std::setw(17) << solution[1]
                << std::setw(17) << solution[2]
                << std::setw(11) << solution[3]
                << std::setw(9) << neu[0]
                << std::setw(9) << neu[1]
                << std::setw(9) << neu[2]
                << std::setw(9) << valid_sats_count
                << std::endl;

            solved_positions.push_back({ solution[0], solution[1], solution[2] });
        }
    }

    calculateRMS(solved_positions, station.approx_pos);
}

void SPPSolver::calculateRMS(const std::vector<Eigen::Vector3d>& positions, const Eigen::Vector3d& ref_pos) {
    if (positions.empty()) {
        std::cout << "\nNo positions to calculate RMS." << std::endl;
        return;
    }

    double sum_dN2 = 0.0, sum_dE2 = 0.0, sum_dU2 = 0.0;
    int n = positions.size();

    Eigen::Vector3d sum_pos = Eigen::Vector3d::Zero();
    for (const auto& pos : positions) {
        sum_pos += pos;
    }
    Eigen::Vector3d avg_pos = sum_pos / n;

    for (const auto& pos : positions) {
        Eigen::Vector3d neu_error = xyz2neu(pos, avg_pos);
        double dN = neu_error[0];
        double dE = neu_error[1];
        double dU = neu_error[2];

        sum_dN2 += dN * dN;
        sum_dE2 += dE * dE;
        sum_dU2 += dU * dU;
    }

    double rms_N = sqrt(sum_dN2 / n);
    double rms_E = sqrt(sum_dE2 / n);
    double rms_U = sqrt(sum_dU2 / n);
    double rms_H = sqrt(rms_N * rms_N + rms_E * rms_E);
    double rms_3D = sqrt(rms_N * rms_N + rms_E * rms_E + rms_U * rms_U);

    std::cout << "\n------------------- Accuracy Assessment (RMS) -------------------" << std::endl;
    std::cout << "Number of Epochs Processed: " << n << std::endl;
    std::cout << "Average Position (XYZ)    : " << avg_pos.transpose() << " (Used as Reference)" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "RMS-North (m)             : " << rms_N << std::endl;
    std::cout << "RMS-East (m)              : " << rms_E << std::endl;
    std::cout << "RMS-Up (m)                : " << rms_U << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << "RMS-Horizontal (2D) (m)   : " << rms_H << std::endl;
    std::cout << "RMS-Total (3D) (m)        : " << rms_3D << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
}

// SPPSolver::processEpoch 方法
bool SPPSolver::processEpoch(double time,
    const std::map<int, Observation>& epoch_obs,
    const Station& station,
    const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris,
    const std::vector<double>& ion_params,
    Eigen::Vector4d& solution,
    int& valid_sats_count) {

    Eigen::Vector3d rec_pos = station.approx_pos;
    double rec_clk_err = 0.0;

    for (int iter = 0; iter < 10; ++iter) {
        std::vector<int> valid_sats;
        std::vector<double> l_vec;
        std::vector<Eigen::Vector4d> A_rows;

        Eigen::Vector3d rec_blh = xyz2blh(rec_pos);

        for (const auto& obs_pair : epoch_obs) {
            int prn = obs_pair.first;
            const Observation& obs = obs_pair.second;

            const GPSEphemeris* eph = selectEphemeris(time, prn, gps_ephemeris);
            if (!eph) continue;

            SatResult sat_res;
            if (!calculateSatPosition(time, obs.C1C, *eph, sat_res)) continue;

            double az, el;
            calculateAzimuthElevation(rec_pos, sat_res.pos, az, el);
            if (el * 180.0 / PI < station.cutoff_angle) continue;

            // 【THE FINAL FIX: TGD CORRECTION】
            double tgd_correction_meters = eph->tgd * C_LIGHT;
            double corrected_pseudorange = obs.C1C - tgd_correction_meters;

            double d_ion = klobucharModel(rec_blh, az, el, time, ion_params);
            double d_trop = saastamoinenModel(el, rec_blh[2]);

            double geometric_dist = (sat_res.pos - rec_pos).norm();

            double l = corrected_pseudorange - (geometric_dist + rec_clk_err - C_LIGHT * sat_res.clock_error + d_ion + d_trop);
            l_vec.push_back(l);

            Eigen::Vector3d e_vec = (rec_pos - sat_res.pos).normalized();
            A_rows.push_back({ e_vec[0], e_vec[1], e_vec[2], 1.0 });

            valid_sats.push_back(prn);
        }

        if (valid_sats.size() < 4) {
            return false;
        }

        valid_sats_count = valid_sats.size();

        Eigen::MatrixXd A(valid_sats.size(), 4);
        Eigen::VectorXd L(valid_sats.size());
        for (size_t i = 0; i < valid_sats.size(); ++i) {
            A.row(i) = A_rows[i];
            L(i) = l_vec[i];
        }

        Eigen::Vector4d dx = (A.transpose() * A).inverse() * A.transpose() * L;

        rec_pos[0] += dx[0];
        rec_pos[1] += dx[1];
        rec_pos[2] += dx[2];
        rec_clk_err += dx[3];

        if (dx.head(3).norm() < 1e-4) {
            solution = { rec_pos[0], rec_pos[1], rec_pos[2], rec_clk_err };
            return true;
        }
    }
    return false;
}

// ... calculateSatPosition, calculateAzimuthElevation, selectEphemeris ...
// (These functions are correct and remain unchanged)
// positioning.cpp

// ... 其他所有函数，包括 solve, calculateRMS, processEpoch 等，都保持我们上一个版本即可 ...

// 只需要替换这一个函数
bool SPPSolver::calculateSatPosition(double obs_time, double pseudorange, const GPSEphemeris& eph, SatResult& result) {
    // 1. 计算信号发射时刻 (用于计算卫星几何位置)
    double tx_time = obs_time - pseudorange / C_LIGHT;

    // 2. 计算卫星位置 (这部分逻辑是正确的)
    double A = eph.sqrtA * eph.sqrtA;
    double n0 = sqrt(GM_EARTH / (A * A * A));
    double n = n0 + eph.delta_n;

    double tk = tx_time - eph.toe;
    if (tk > 302400) tk -= 604800;
    if (tk < -302400) tk += 604800;

    double Mk = eph.m0 + n * tk;

    double Ek = Mk;
    for (int i = 0; i < 10; ++i) {
        Ek = Mk + eph.e * sin(Ek);
    }

    double vk = atan2(sqrt(1.0 - eph.e * eph.e) * sin(Ek), cos(Ek) - eph.e);
    double phik = vk + eph.w;

    double duk = eph.cus * sin(2 * phik) + eph.cuc * cos(2 * phik);
    double drk = eph.crs * sin(2 * phik) + eph.crc * cos(2 * phik);
    double dik = eph.cis * sin(2 * phik) + eph.cic * cos(2 * phik);

    double uk = phik + duk;
    double rk = A * (1.0 - eph.e * cos(Ek)) + drk;
    double ik = eph.i0 + eph.idot * tk + dik;

    double xk_prime = rk * cos(uk);
    double yk_prime = rk * sin(uk);

    double omegak = eph.omega0 + (eph.omega_dot - OMEGA_EARTH) * tk - OMEGA_EARTH * eph.toe;

    result.pos[0] = xk_prime * cos(omegak) - yk_prime * cos(ik) * sin(omegak);
    result.pos[1] = xk_prime * sin(omegak) + yk_prime * cos(ik) * cos(omegak);
    result.pos[2] = yk_prime * sin(ik);

    // 3. 【最终修正】计算卫星钟差 (必须使用信号接收时刻 obs_time 作为参考)
    double dt = obs_time - eph.toc; // <--- 使用 obs_time，而不是 tx_time
    // 同样需要处理周翻转
    if (dt > 302400) dt -= 604800;
    if (dt < -302400) dt += 604800;

    // 多项式部分
    result.clock_error = eph.af0 + eph.af1 * dt + eph.af2 * dt * dt;

    // 相对论效应改正 (F = -2*sqrt(GM)/c^2)
    double rel_effect = -4.442807633e-10 * eph.e * eph.sqrtA * sin(Ek);
    result.clock_error += rel_effect;

    return true;
}

void SPPSolver::calculateAzimuthElevation(const Eigen::Vector3d& rec_pos, const Eigen::Vector3d& sat_pos, double& az, double& el) {
    Eigen::Vector3d d_xyz = sat_pos - rec_pos;
    Eigen::Vector3d neu = xyz2neu(sat_pos, rec_pos);
    el = asin(neu[2] / d_xyz.norm());
    az = atan2(neu[1], neu[0]);
    if (az < 0) az += 2 * PI;
}

const GPSEphemeris* SPPSolver::selectEphemeris(double time, int prn, const std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris) {
    if (gps_ephemeris.find(prn) == gps_ephemeris.end()) return nullptr;
    const GPSEphemeris* best_eph = nullptr;
    double min_diff = 1e9;
    for (const auto& eph : gps_ephemeris.at(prn)) {
        double diff = abs(time - eph.toe);
        if (diff < min_diff) {
            min_diff = diff;
            best_eph = &eph;
        }
    }
    return best_eph;
}