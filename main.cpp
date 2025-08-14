#include "rinex.h"
#include "positioning.h"
#include <iostream>

int main() {
    // 1. �����ļ�·�� 
    std::string obs_filepath = "E:\\GNSS_Date\\wuhn0020.24o";
    std::string nav_filepath = "E:\\GNSS_Date\\brdm0020.24p";

    // 2. ��ʼ��
    RinexParser parser;
    SPPSolver solver;

    Station station;
    station.cutoff_angle = 15.0; // ���ý�ֹ�߶Ƚ�

    std::map<int, std::vector<GPSEphemeris>> gps_ephemeris;
    std::vector<double> ion_params;
    std::map<double, std::map<int, Observation>> obs_data;

    // 3. ��ȡ�ļ� (��Ӧ����ͼ "��ȡ�����ļ�" �� "��ȡ�۲��ļ�")
    std::cout << "Reading navigation file..." << std::endl;
    if (!parser.readNavFile(nav_filepath, gps_ephemeris, ion_params)) {
        return -1;
    }
    std::cout << "Navigation file read successfully." << std::endl;

    std::cout << "Reading observation file..." << std::endl;
    if (!parser.readObsFile(obs_filepath, station, obs_data)) {
        return -1;
    }
    std::cout << "Observation file read successfully." << std::endl;
    std::cout << "Station approximate position (X,Y,Z): "
        << station.approx_pos.transpose() << std::endl;

    // 4. ִ��SPP���� (��Ӧ����ͼ�ĺ���ѭ������)
    std::cout << "\nStarting SPP calculation..." << std::endl;
    solver.solve(station, obs_data, gps_ephemeris, ion_params);

    std::cout << "\nSPP processing finished." << std::endl;

    system("pause");
    return 0;
}