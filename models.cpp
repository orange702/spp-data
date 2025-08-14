#include "models.h"
#include <cmath>

double klobucharModel(const Eigen::Vector3d& rec_pos_blh, double az, double el, double gpst, const std::vector<double>& ion_params) {
    if (ion_params.size() < 8) return 0.0;
    double lat = rec_pos_blh[0];
    double lon = rec_pos_blh[1];

    double psi = 0.0137 / (el / PI + 0.11) - 0.022;
    double phi_i = lat / PI + psi * cos(az);
    if (phi_i > 0.416) phi_i = 0.416;
    if (phi_i < -0.416) phi_i = -0.416;

    double lambda_i = lon / PI + (psi * sin(az)) / cos(phi_i * PI);
    double phi_m = phi_i + 0.064 * cos((lambda_i - 1.617) * PI);

    double t = fmod(4.32e4 * lambda_i + gpst, 86400.0);
    if (t > 86400.0) t -= 86400.0;
    if (t < 0) t += 86400.0;

    double PER = ion_params[4] + ion_params[5] * phi_m + ion_params[6] * pow(phi_m, 2) + ion_params[7] * pow(phi_m, 3);
    if (PER < 72000.0) PER = 72000.0;

    double AMP = ion_params[0] + ion_params[1] * phi_m + ion_params[2] * pow(phi_m, 2) + ion_params[3] * pow(phi_m, 3);
    if (AMP < 0) AMP = 0.0;

    double x = 2.0 * PI * (t - 50400.0) / PER;
    double F = 1.0 + 16.0 * pow(0.53 - el / PI, 3.0);
    double d_ion = 0.0;
    if (abs(x) < 1.57) {
        d_ion = F * (5e-9 + AMP * (1.0 - pow(x, 2) / 2.0 + pow(x, 4) / 24.0));
    }
    else {
        d_ion = F * 5e-9;
    }
    return C_LIGHT * d_ion;
}

double saastamoinenModel(double el, double rec_height) {
    if (el <= 0) return 0.0;
    double P = 1013.25 * pow(1 - 2.2557e-5 * rec_height, 5.2559);
    double T = 15.0 - 0.0065 * rec_height + 273.15;
    double e = 0.5 * exp(24.37 - 6142.0 / (273.0 + T - 273.15));

    double z = PI / 2.0 - el;
    double d_ztd = 0.002277 / cos(z) * (P - (0.0005 + 0.002277 * (1255.0 / T + 0.05)) * e);
    return d_ztd;
}