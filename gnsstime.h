#pragma once

// 将公历时间转换为GPS周和周内秒
void ymdhmsToGPS(int year, int month, int day, int hour, int min, double sec, int& week, double& sow);

// 将GPS周和周内秒转换为公历时间
void gpsToYmdhms(int week, double sow, int& year, int& month, int& day, int& hour, int& min, double& sec);