#pragma once

#include "global.h"
#include <iostream>
#include <fstream>
#include <sstream>

class RinexParser {
public:
    // ��ȡ�����ļ�
    bool readNavFile(const std::string& filepath, std::map<int, std::vector<GPSEphemeris>>& gps_ephemeris, std::vector<double>& ion_params);

    // ��ȡ�۲��ļ�
    bool readObsFile(const std::string& filepath, Station& station, std::map<double, std::map<int, Observation>>& obs_data);

private:
    // ����RINEX 3.xx �����ļ�ͷ��
    void parseNavHeader(std::ifstream& file, std::vector<double>& ion_params);

    // ����RINEX 3.xx �۲��ļ�ͷ��
    void parseObsHeader(std::ifstream& file, Station& station);
};