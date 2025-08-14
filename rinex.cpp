#include "rinex.h"
#include "gnsstime.h"

bool RinexParser::readNavFile(const std::string& filepath, std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris, std::vector<double>& ion_params) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open navigation file: " << filepath << std::endl;
        return false;
    }

    std::string line;
    // 读取文件头
    while (getline(file, line) && line.find("END OF HEADER") == std::string::npos) {
        if (line.find("IONOSPHERIC CORR") != std::string::npos && line.substr(0, 4) == "GPSA") {
            for (int i = 0; i < 4; ++i) ion_params.push_back(std::stod(line.substr(5 + i * 12, 12)));
        }
        if (line.find("IONOSPHERIC CORR") != std::string::npos && line.substr(0, 4) == "GPSB") {
            for (int i = 0; i < 4; ++i) ion_params.push_back(std::stod(line.substr(5 + i * 12, 12)));
        }
    }

    // 读取数据体
    while (getline(file, line)) {
        if (line[0] != 'G') continue; // 只处理GPS卫星

        GPSEphemeris eph;
        eph.prn = std::stoi(line.substr(1, 2));

        // 时间
        int year = std::stoi(line.substr(4, 4));
        int month = std::stoi(line.substr(9, 2));
        int day = std::stoi(line.substr(12, 2));
        int hour = std::stoi(line.substr(15, 2));
        int min = std::stoi(line.substr(18, 2));
        double sec = std::stod(line.substr(21, 2));
        int week;
        ymdhmsToGPS(year, month, day, hour, min, sec, week, eph.toc);
        eph.week = week;

        // 解析参数
        eph.af0 = std::stod(line.substr(23, 19));
        eph.af1 = std::stod(line.substr(42, 19));
        eph.af2 = std::stod(line.substr(61, 19));

        getline(file, line); eph.iode = std::stod(line.substr(4, 19)); eph.crs = std::stod(line.substr(23, 19)); eph.delta_n = std::stod(line.substr(42, 19)); eph.m0 = std::stod(line.substr(61, 19));
        getline(file, line); eph.cuc = std::stod(line.substr(4, 19)); eph.e = std::stod(line.substr(23, 19)); eph.cus = std::stod(line.substr(42, 19)); eph.sqrtA = std::stod(line.substr(61, 19));
        getline(file, line); eph.toe = std::stod(line.substr(4, 19)); eph.cic = std::stod(line.substr(23, 19)); eph.omega0 = std::stod(line.substr(42, 19)); eph.cis = std::stod(line.substr(61, 19));
        getline(file, line); eph.i0 = std::stod(line.substr(4, 19)); eph.crc = std::stod(line.substr(23, 19)); eph.w = std::stod(line.substr(42, 19)); eph.omega_dot = std::stod(line.substr(61, 19));
        getline(file, line); eph.idot = std::stod(line.substr(4, 19)); // week is already parsed from line 6.
        getline(file, line); eph.health = (std::stod(line.substr(23, 19)) == 0.0); eph.tgd = std::stod(line.substr(42, 19));
        getline(file, line); // 最后一行不含关键信息，跳过

        gps_ephemeris[eph.prn].push_back(eph);
    }
    return true;
}

bool RinexParser::readObsFile(const std::string& filepath, Station& station, std::map<double, std::map<int, Observation>>& obs_data) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open observation file: " << filepath << std::endl;
        return false;
    }

    std::string line;
    // 读取文件头
    while (getline(file, line) && line.find("END OF HEADER") == std::string::npos) {
        if (line.find("APPROX POSITION XYZ") != std::string::npos) {
            station.approx_pos[0] = std::stod(line.substr(0, 14));
            station.approx_pos[1] = std::stod(line.substr(14, 14));
            station.approx_pos[2] = std::stod(line.substr(28, 14));
        }
    }

    // 读取数据体
    while (getline(file, line)) {
        if (line[0] != '>') continue; // 历元开始标志

        int year = std::stoi(line.substr(2, 4));
        int month = std::stoi(line.substr(7, 2));
        int day = std::stoi(line.substr(10, 2));
        int hour = std::stoi(line.substr(13, 2));
        int min = std::stoi(line.substr(16, 2));
        double sec = std::stod(line.substr(19, 10));
        int num_sats = std::stoi(line.substr(32, 3));

        int week;
        double sow;
        ymdhmsToGPS(year, month, day, hour, min, sec, week, sow);

        for (int i = 0; i < num_sats; ++i) {
            getline(file, line);
            if (line[0] != 'G') continue; // 只处理GPS

            int prn = std::stoi(line.substr(1, 2));
            try {
                double c1c = std::stod(line.substr(3, 14));
                obs_data[sow][prn].time = sow;
                obs_data[sow][prn].C1C = c1c;
            }
            catch (const std::invalid_argument& e) {
                // 忽略没有C1C观测值的行
            }
        }
    }
    return true;
}