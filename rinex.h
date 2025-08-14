#pragma once

#include "global.h"
#include <iostream>
#include <fstream>
#include <sstream>

class RinexParser {
public:
    // 读取导航文件
    bool readNavFile(const std::string& filepath, std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris, std::vector<double>& ion_params);

    // 读取观测文件
    bool readObsFile(const std::string& filepath, Station& station, std::map<double, std::map<int, Observation>>& obs_data);

private:
    // 解析RINEX 3.xx 导航文件头部
    void parseNavHeader(std::ifstream& file, std::vector<double>& ion_params);

    // 解析RINEX 3.xx 观测文件头部
    void parseObsHeader(std::ifstream& file, Station& station);
};